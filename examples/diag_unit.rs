use mfs::design::FilterDesign;
use num_complex::Complex64;
fn pv(c:&[Complex64],s:Complex64)->Complex64{c.iter().rev().fold(Complex64::new(0.0,0.0),|a,x|a*s+x)}
fn alt(c:&[Complex64])->Vec<Complex64>{c.iter().enumerate().map(|(k,x)|x.conj()*Complex64::new((-1.0f64).powi(k as i32),0.0)).collect()}
fn conv(a:&[Complex64],b:&[Complex64])->Vec<Complex64>{let mut o=vec![Complex64::new(0.0,0.0);a.len()+b.len()-1];for i in 0..a.len(){for j in 0..b.len(){o[i+j]+=a[i]*b[j];}}o}
fn main(){
    let d=FilterDesign::bandpass(2,22.0,1.95e9,200e6).synthesize().unwrap();
    let p=d.polynomials();
    let e=&p.e.coefficients;let f=&p.f.coefficients;let pp=&p.p.coefficients;
    let diff=conv(e,&alt(e));let ff=conv(f,&alt(f));let ppp=conv(pp,&alt(pp));
    println!("E=", );
    for (i,c) in e.iter().enumerate(){print!("e{i}=({:.4},{:.4}) ",c.re,c.im);}println!();
    for (i,c) in f.iter().enumerate(){print!("f{i}=({:.4},{:.4}) ",c.re,c.im);}println!();
    for (i,c) in pp.iter().enumerate(){print!("p{i}=({:.4},{:.4}) ",c.re,c.im);}println!();
    println!("EE~ - FF~ - PP~ coeffs:");
    for i in 0..diff.len(){let r=diff[i]-ff.get(i).copied().unwrap_or_default()-ppp.get(i).copied().unwrap_or_default();print!("  s^{i}: ({:.4},{:.4})",r.re,r.im);}println!();
    for w in [-1.0,-0.5,0.0,0.5,1.0]{
        let s=Complex64::new(0.0,w);
        let ev=pv(e,s);let fv=pv(f,s);let pvv=pv(pp,s);
        println!(" w={w:+.1} |S11|^2+|S21|^2 = {:.10}",fv.norm_sqr()/ev.norm_sqr()+pvv.norm_sqr()/ev.norm_sqr());
    }
}
