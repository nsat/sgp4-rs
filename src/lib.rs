use chrono::prelude::*;
use chrono::DateTime;

use thiserror::Error;

mod sgp4_sys;

#[derive(Debug, Error)]
pub enum Error {
    #[error("TLE was malformed: {0}")]
    MalformedTwoLineElement(String),
    #[error("Error in SGP4 propagator")]
    PropagationError,
    #[error("{0}")]
    UnknownError(String),
}

type Result<T> = std::result::Result<T, Error>;

pub struct StateVector {
    pub position: [f64; 3],
    pub velocity: [f64; 3],
}

const TLE_LINE_LENGTH: usize = 69;

pub struct TwoLineElement {
    elements: sgp4_sys::OrbitalElementSet,
}

impl TwoLineElement {
    /// Create a validated TwoLineElement from a string.
    pub fn new(line1: &str, line2: &str) -> Result<TwoLineElement> {
        if line1.len() != TLE_LINE_LENGTH || line2.len() != TLE_LINE_LENGTH {
            return Err(Error::MalformedTwoLineElement("Line is the wrong length".to_string()));
        }

        let elements = sgp4_sys::to_orbital_elements(
            line1,
            line2,
            sgp4_sys::RunType::Verification,
            sgp4_sys::OperationMode::Improved,
            sgp4_sys::GravitationalConstant::Wgs84,
        )
        .map_err(|e| Error::MalformedTwoLineElement(format!("{:?}", e)))?;

        Ok(TwoLineElement { elements })
    }

    /// Get the epoch of a TwoLineElement.
    pub fn epoch(&self) -> Result<DateTime<Utc>> {
        Ok(self.elements.epoch())
    }

    pub fn propagate_to(&self, t: DateTime<Utc>) -> Result<StateVector> {
        let tle_epoch = self.elements.epoch();
        let min_since_epoch = (t - tle_epoch).num_days() as f64;

        let (r, v) = sgp4_sys::run_sgp4(
            self.elements,
            sgp4_sys::GravitationalConstant::Wgs84,
            min_since_epoch,
        )
        .map_err(|_e| Error::PropagationError)?;

        Ok(StateVector {
            position: r.to_owned(),
            velocity: v.to_owned(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use chrono::Duration;

    #[test]
    fn test_simple_propagation() -> Result<()> {
        let line1 = "1 25544U 98067A   20148.21301450  .00001715  00000-0  38778-4 0  9992";
        let line2 = "2 25544  51.6435  92.2789 0002570 358.0648 144.9972 15.49396855228767";

        let tle = TwoLineElement::new(line1, line2)?;
        let epoch = tle.epoch()?;

        let _s1 = tle.propagate_to(epoch);
        let _s2 = tle.propagate_to(epoch + Duration::hours(1));

        Ok(())
    }
}
