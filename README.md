# sgp4-rs

[![CI Status](https://github.com/nsat/sgp4-rs/workflows/Rust/badge.svg)](https://github.com/nsat/sgp4-rs/actions)
[![Docs](https://docs.rs/sgp4-rs/badge.svg)](https://docs.rs/sgp4-rs/)
[![Crate](https://img.shields.io/crates/v/sgp4-rs)](https://crates.io/crates/sgp4-rs)

This crate implements a wrapper around the C++ implementation of the SGP-4 orbital propagator as
provided in the "Revisiting Spacetrack Report #3" paper
([link](https://celestrak.com/publications/AIAA/2006-6753/)). It provides high level bindings to the
propagator library with a more modern interface.

Our approach separates the low-level `unsafe` bindings into the `sgp4_sys` module, while safe
functions are exported through the library's root module. Because the underlying SGP4 implementation
is thread-safe, this crate can be used in multithreaded environments and with async/await code.

We have not created bindings to every function in the library, especially as some of them are
duplicative of Rust standard library functions. The core propagator functionality is exposed, and
allows predicting an orbiting body's state vector at a given time from two line element data.

## Building

`sgp4` builds cleanly on the stable Rust channel, but does require a local C++ compiler to be
present in order to build the wrapped SGP4 library.

## Experimental Features

The `tlegen` feature adds basic support for creating custom TLEs from a set of orbital elements.
This feature is subject to several important caveats, so it is not enabled by default. See the
`ClassicalOrbitalElements` documentation for details.

## Related

The [`sgp4` crate](https://github.com/neuromorphicsystems/sgp4) is a pure Rust reimplementation of
the Vallado library which this crate wraps.
