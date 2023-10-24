# Building checks on non-standard configurations


### Install docker 

From terminal, run:
```
  $ sudo apt-get install -y docker.io

  $ sudo usermod -aG docker $USER
```
Finally, reboot your PC or run, from terminal, the following command:
  
```
  $ newgrp docker
```

### Running building checks
Different shell scripts are available to check the development version of `fdaPDE` package on non-standard configurations. Assuming your working directory is `fdaPDE/tests/building_checks/`, run, from terminal, one of the following commands.

  - Building check on Fedora using R-devel and gcc as compiler, run from terminal:
  ```
    $ ./fedora_gcc_devel.sh
  ```
  
  - Building check on Fedora using R-devel and clang as compiler, run from terminal:
  ```
    $ ./fedora_clang_devel.sh
  ```
  
  - Building check on Debian using R-patched and gcc as compiler, run from terminal:
  ```
    $ ./debian_gcc_patched.sh
  ```

Finally, check the directory `fdaPDE/fdaPDE.Rcheck/`.
  
