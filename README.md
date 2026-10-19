# compressibleTools
OpenFOAM library for compressible flow calculations - boundary conditions, linear solvers (BiCGStab from extend-project), function objects

## Building

This branch is laid out as an OpenFOAM (openfoam.com) plugin.

- Standalone build into user locations (`$FOAM_USER_LIBBIN`):
  ```bash
  ./Allwmake
  ```
- As a plugin: place (or add as a git submodule) this repository in
  `$WM_PROJECT_DIR/plugins/libcompressibleTools` and run `./Allwmake-plugins` from
  `$WM_PROJECT_DIR`. It is then installed into `$FOAM_LIBBIN`.
- Custom location:
  `./Allwmake -prefix=user|group|openfoam|/some/path`
