# Current release process for Rust crate sgp4-rs

- Get any outstanding PRs you want in the release merged to master.
- Ensure you have an account on https://crates.io/ linked to your
  github account.
- Make sure your `crates.io` account has a verified email. 
- Use `cargo login` to log in to your crates.io account locally.
- Checkout/pull the repo master branch locally
- Check you are one of the owners of the crate, using `cargo owner
  --list sgp4-rs`. If not, persuade one of the other owners to add
  you, and make sure you accept that ownership by clicking in the email.
- Choose whether you are doing a major, minor or patch release.
- Assuming patch, within the root dir of the repo, do `cargo release
  patch --execute`.
- This should bump the version in the `Cargo.toml` file, push that
  update, tag the repo appropriately and publish the crate to
  `crates.io`
- Note, this does not create a GitHub release.
