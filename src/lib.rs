use chrono::prelude::*;
use chrono::{DateTime, Duration};

use thiserror::Error;

mod sgp4_sys;

#[derive(Debug, Error)]
pub enum Error {
    #[error("TLE was malformed")]
    MalformedTwoLineElement,
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

const TLE_LINE_LENGTH: usize = 70;

pub struct TwoLineElement {
    line1: String,
    line2: String,
}

impl TwoLineElement {
    pub fn new(line1: &str, line2: &str) -> Result<TwoLineElement> {
        if line1.len() != TLE_LINE_LENGTH || line2.len() != TLE_LINE_LENGTH {
            return Err(Error::MalformedTwoLineElement);
        }
        Ok(TwoLineElement { line1: line1.to_owned(), line2: line2.to_owned() })
    }

    pub fn epoch(&self) -> Result<DateTime<Utc>> {
        let sat_state = self
            .to_orbital_elements()
            .map_err(|_e| Error::MalformedTwoLineElement)?;

        Ok(sat_state.epoch())
    }

    pub fn propagate_to(&self, dt: DateTime<Utc>) -> Result<StateVector> {
        let sat_state = self
            .to_orbital_elements()
            .map_err(|_e| Error::MalformedTwoLineElement)?;

        let tle_epoch = sat_state.epoch();
        let min_since_epoch = (dt - tle_epoch).num_days() as f64;

        let (r, v) = sgp4_sys::run_sgp4(
            sat_state,
            sgp4_sys::GravitationalConstant::Wgs84,
            min_since_epoch,
        )
        .map_err(|_e| Error::PropagationError)?;

        Ok(StateVector {
            position: r.to_owned(),
            velocity: v.to_owned(),
        })
    }

    fn to_orbital_elements(&self) -> Result<sgp4_sys::OrbitalElementSet> {
        sgp4_sys::to_orbital_elements(
            &self.line1,
            &self.line2,
            sgp4_sys::RunType::Verification,
            sgp4_sys::OperationMode::Improved,
            sgp4_sys::GravitationalConstant::Wgs84,
        )
        .map_err(|e| Error::UnknownError(format!("{:?}", e)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_propagation() -> Result<()> {
        let line1 = "1 25544U 98067A   20148.21301450  .00001715  00000-0  38778-4 0  9992";
        let line2 = "2 25544  51.6435  92.2789 0002570 358.0648 144.9972 15.49396855228767";

        let tle = TwoLineElement::new(line1, line2)?;
        let epoch = tle.epoch()?;

        let s1 = tle.propagate_to(epoch);
        let s2 = tle.propagate_to(epoch + Duration::hours(1));

        Ok(())
    }
}
