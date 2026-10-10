mod builder;
mod multi_band;
mod types;

pub use builder::FilterSpecBuilder;
pub use multi_band::{BandSpec, MultiBandSpec, RejectionRegion, RejectionSpec};
pub use types::{
    FilterSpec, MAX_ORDER, NormalizedFilterSpec, OutOfBandAttenuationSpec,
    OutOfBandAttenuationWindow, TransmissionZero,
};
