#!/usr/bin/env python3
"""
Extract filter design data from PDF papers for MFS test fixtures.

Usage:
    python extract_from_pdf.py paper.pdf --output fixtures.json
    python extract_from_pdf.py paper.pdf --pages 5-8 --output fixtures.json
    python extract_from_pdf.py paper.pdf --interactive

Requirements:
    pip install pymupdf tabula-py pandas

This script extracts:
- Filter specifications (order, return loss, transmission zeros)
- Polynomial coefficients (E, F, P)
- Coupling matrix values
- S-parameter reference data points
"""

import argparse
import json
import re
import sys
from pathlib import Path

try:
    import fitz  # PyMuPDF
except ImportError:
    fitz = None
    print("Warning: pymupdf not installed. Install with: pip install pymupdf")

try:
    import tabula
except ImportError:
    tabula = None
    print("Warning: tabula-py not installed. Install with: pip install tabula-py")


# ─── Pattern matchers for common filter paper formats ─────────────────────────

# Matches: "order 4", "N = 6", "4th-order", "6-pole"
ORDER_PATTERNS = [
    r'(?:order|N)\s*[=:]\s*(\d+)',
    r'(\d+)(?:th|rd|nd|st)[\s-]*order',
    r'(\d+)[\s-]*pole',
]

# Matches: "return loss 20 dB", "RL = 22 dB", "20-dB return loss"
RETURN_LOSS_PATTERNS = [
    r'(?:return\s*loss|RL)\s*[=:]\s*([\d.]+)\s*dB',
    r'([\d.]+)[\s-]*dB\s*return\s*loss',
]

# Matches transmission zeros: "±1.5", "1.3217j", "jΩ = ±1.5"
TZ_PATTERNS = [
    r'[±]\s*([\d.]+)',
    r'(?:transmission\s*zero|TZ).*?([\d.]+)',
    r'j[Ωω]\s*=\s*[±]?\s*([\d.]+)',
]

# Matches epsilon values: "ε = 1.1548", "εR = 1.0"
EPSILON_PATTERNS = [
    r'[εe](?:psilon)?\s*[=:]\s*([\d.]+)',
    r'[εe](?:psilon)?_?[Rr]\s*[=:]\s*([\d.]+)',
]

# Matches matrix elements: "M12 = 0.8832", "m_{12} = 0.8832"
MATRIX_ELEMENT_PATTERN = r'[Mm]_?\{?(\d)(\d)\}?\s*[=:]\s*([+-]?[\d.]+)'


def extract_text_from_pdf(pdf_path, pages=None):
    """Extract text from PDF using PyMuPDF."""
    if fitz is None:
        raise ImportError("pymupdf is required: pip install pymupdf")
    
    doc = fitz.open(pdf_path)
    texts = []
    
    if pages:
        page_range = parse_page_range(pages, len(doc))
    else:
        page_range = range(len(doc))
    
    for page_num in page_range:
        page = doc[page_num]
        texts.append({
            'page': page_num + 1,
            'text': page.get_text(),
        })
    
    doc.close()
    return texts


def extract_tables_from_pdf(pdf_path, pages=None):
    """Extract tables from PDF using tabula-py."""
    if tabula is None:
        raise ImportError("tabula-py is required: pip install tabula-py")
    
    kwargs = {'pages': pages or 'all', 'multiple_tables': True}
    try:
        tables = tabula.read_pdf(pdf_path, **kwargs)
        return tables
    except Exception as e:
        print(f"Warning: tabula extraction failed: {e}")
        return []


def parse_page_range(pages_str, total_pages):
    """Parse page range string like '5-8' or '1,3,5-7'."""
    result = []
    for part in pages_str.split(','):
        if '-' in part:
            start, end = part.split('-')
            result.extend(range(int(start) - 1, min(int(end), total_pages)))
        else:
            result.append(int(part) - 1)
    return result


def find_filter_specs(text):
    """Extract filter specifications from text."""
    specs = []
    
    # Find order
    order = None
    for pattern in ORDER_PATTERNS:
        match = re.search(pattern, text, re.IGNORECASE)
        if match:
            order = int(match.group(1))
            break
    
    # Find return loss
    return_loss = None
    for pattern in RETURN_LOSS_PATTERNS:
        match = re.search(pattern, text, re.IGNORECASE)
        if match:
            return_loss = float(match.group(1))
            break
    
    # Find transmission zeros
    zeros = []
    for pattern in TZ_PATTERNS:
        matches = re.findall(pattern, text, re.IGNORECASE)
        for m in matches:
            val = float(m)
            if 0.5 < val < 10.0:  # Reasonable range for normalized zeros
                if val not in zeros:
                    zeros.append(val)
    
    # Find epsilon
    epsilon = None
    epsilon_r = None
    for pattern in EPSILON_PATTERNS:
        matches = re.findall(pattern, text)
        if matches:
            if epsilon is None:
                epsilon = float(matches[0])
            if len(matches) > 1 and epsilon_r is None:
                epsilon_r = float(matches[1])
    
    if order or return_loss or zeros:
        specs.append({
            'order': order,
            'return_loss_db': return_loss,
            'transmission_zeros': zeros,
            'epsilon': epsilon,
            'epsilon_r': epsilon_r,
        })
    
    return specs


def find_matrix_elements(text):
    """Extract coupling matrix elements from text."""
    elements = {}
    for match in re.finditer(MATRIX_ELEMENT_PATTERN, text):
        row, col, value = int(match.group(1)), int(match.group(2)), float(match.group(3))
        elements[(row, col)] = value
    return elements


def find_polynomial_coefficients(text):
    """Extract polynomial coefficients from text (E, F, P polynomials)."""
    polys = {}
    
    # Look for coefficient lists after E(s), F(s), P(s) labels
    for poly_name in ['E', 'F', 'P']:
        # Pattern: E(s) = ... or E: [...]
        pattern = rf'{poly_name}\s*(?:\(s\))?\s*[=:]\s*\[?([\d.,\s+-]+)\]?'
        match = re.search(pattern, text)
        if match:
            coeff_str = match.group(1)
            try:
                coeffs = [float(x.strip()) for x in re.split(r'[,\s]+', coeff_str) if x.strip()]
                if coeffs:
                    polys[poly_name] = coeffs
            except ValueError:
                pass
    
    return polys


def build_fixture(spec, matrix_elements, polynomials, source_info):
    """Build a test fixture JSON object from extracted data."""
    fixture = {
        'case_id': f"extracted_{spec.get('order', 'N')}_{len(spec.get('transmission_zeros', []))}TZ",
        'source': source_info,
        'specification': {
            'filter_order': spec.get('order'),
            'return_loss': {'value': spec.get('return_loss_db'), 'unit': 'dB'},
        },
    }
    
    if spec.get('transmission_zeros'):
        fixture['specification']['normalized_transmission_zeros'] = [
            {'value': {'re': 0.0, 'im': z}, 'unit': 'rad/s', 'type': 'pure_imag', 'domain': 's'}
            for z in spec['transmission_zeros']
        ]
    
    if polynomials:
        fixture['mathematical_model'] = {
            'domain': 's',
            'polynomial_coefficients': polynomials,
        }
        if spec.get('epsilon') or spec.get('epsilon_r'):
            fixture['mathematical_model']['singularities'] = {}
            if spec.get('epsilon'):
                fixture['mathematical_model']['singularities']['epsilon'] = spec['epsilon']
            if spec.get('epsilon_r'):
                fixture['mathematical_model']['singularities']['epsilon_R'] = spec['epsilon_r']
    
    if matrix_elements:
        # Reconstruct matrix from elements
        max_idx = max(max(r, c) for r, c in matrix_elements.keys())
        matrix = [[0.0] * (max_idx + 1) for _ in range(max_idx + 1)]
        for (r, c), v in matrix_elements.items():
            matrix[r][c] = v
            matrix[c][r] = v  # symmetric
        fixture['coupling_matrix'] = {
            'topology': 'folded',
            'data': matrix,
        }
    
    return fixture


def process_pdf(pdf_path, pages=None):
    """Process a PDF and extract all filter design data."""
    print(f"Processing: {pdf_path}")
    
    # Extract text
    page_texts = extract_text_from_pdf(pdf_path, pages)
    full_text = '\n'.join(p['text'] for p in page_texts)
    
    # Extract data
    specs = find_filter_specs(full_text)
    matrix_elements = find_matrix_elements(full_text)
    polynomials = find_polynomial_coefficients(full_text)
    
    # Source info
    source_info = {
        'file': str(pdf_path),
        'pages_scanned': len(page_texts),
    }
    
    # Try to extract title from first page
    if page_texts:
        first_lines = page_texts[0]['text'].split('\n')[:5]
        title_candidate = max(first_lines, key=len) if first_lines else ''
        if len(title_candidate) > 10:
            source_info['title'] = title_candidate.strip()
    
    # Build fixtures
    fixtures = []
    for spec in specs:
        fixture = build_fixture(spec, matrix_elements, polynomials, source_info)
        fixtures.append(fixture)
    
    # Also try table extraction
    tables = extract_tables_from_pdf(pdf_path, pages)
    if tables:
        print(f"  Found {len(tables)} tables")
        for i, table in enumerate(tables):
            print(f"  Table {i+1}: {table.shape[0]} rows x {table.shape[1]} cols")
            # Try to identify coupling matrix tables
            if table.shape[0] == table.shape[1] and 3 <= table.shape[0] <= 12:
                print(f"    → Possible coupling matrix ({table.shape[0]}x{table.shape[0]})")
    
    return fixtures, tables


def interactive_mode(pdf_path):
    """Interactive mode for manual data entry guided by PDF content."""
    print(f"\n{'='*60}")
    print(f"Interactive Filter Data Extraction")
    print(f"PDF: {pdf_path}")
    print(f"{'='*60}\n")
    
    fixture = {
        'case_id': input("Case ID (e.g., Cameron2003_Table1): ").strip(),
        'source': {
            'title': input("Paper title: ").strip(),
            'authors': input("Authors: ").strip(),
            'year': int(input("Year: ").strip() or "2003"),
        },
        'specification': {},
    }
    
    order = input("Filter order: ").strip()
    if order:
        fixture['specification']['filter_order'] = int(order)
    
    rl = input("Return loss (dB): ").strip()
    if rl:
        fixture['specification']['return_loss'] = {'value': float(rl), 'unit': 'dB'}
    
    zeros_str = input("Transmission zeros (comma-separated, e.g., -1.5,1.5): ").strip()
    if zeros_str:
        zeros = [float(z) for z in zeros_str.split(',')]
        fixture['specification']['normalized_transmission_zeros'] = [
            {'value': {'re': 0.0, 'im': z}, 'unit': 'rad/s', 'type': 'pure_imag', 'domain': 's'}
            for z in zeros
        ]
    
    eps = input("Epsilon (or Enter to skip): ").strip()
    if eps:
        fixture['mathematical_model'] = {
            'singularities': {'epsilon': float(eps)}
        }
        eps_r = input("Epsilon_R (or Enter for 1.0): ").strip()
        fixture['mathematical_model']['singularities']['epsilon_R'] = float(eps_r) if eps_r else 1.0
    
    matrix_str = input("\nCoupling matrix (paste rows separated by ';', or Enter to skip):\n").strip()
    if matrix_str:
        rows = matrix_str.split(';')
        matrix = [[float(x) for x in row.split(',')] for row in rows]
        fixture['coupling_matrix'] = {
            'topology': input("Topology (folded/arrow/transversal): ").strip() or 'folded',
            'data': matrix,
        }
    
    return fixture


def main():
    parser = argparse.ArgumentParser(description='Extract filter design data from PDF papers')
    parser.add_argument('pdf', help='Path to PDF file')
    parser.add_argument('--output', '-o', default='extracted_fixtures.json', help='Output JSON file')
    parser.add_argument('--pages', '-p', help='Page range (e.g., "5-8" or "1,3,5-7")')
    parser.add_argument('--interactive', '-i', action='store_true', help='Interactive data entry mode')
    parser.add_argument('--append', '-a', action='store_true', help='Append to existing output file')
    args = parser.parse_args()
    
    pdf_path = Path(args.pdf)
    if not pdf_path.exists():
        print(f"Error: {pdf_path} not found")
        sys.exit(1)
    
    if args.interactive:
        fixture = interactive_mode(pdf_path)
        fixtures = [fixture]
        tables = []
    else:
        fixtures, tables = process_pdf(pdf_path, args.pages)
    
    # Load existing data if appending
    output_path = Path(args.output)
    if args.append and output_path.exists():
        with open(output_path) as f:
            existing = json.load(f)
        if isinstance(existing, dict) and 'case' in existing:
            existing['case'].extend(fixtures)
            output_data = existing
        elif isinstance(existing, list):
            existing.extend(fixtures)
            output_data = existing
        else:
            output_data = {'schema_version': '1.0', 'case': fixtures}
    else:
        output_data = {
            'schema_version': '1.0',
            'extracted_from': str(pdf_path),
            'case': fixtures,
        }
    
    # Write output
    with open(output_path, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"\nExtracted {len(fixtures)} fixture(s) → {output_path}")
    
    if not fixtures:
        print("\nNo filter data automatically detected.")
        print("Try: python extract_from_pdf.py paper.pdf --interactive")
        print("Or manually add data to the JSON file.")


if __name__ == '__main__':
    main()
