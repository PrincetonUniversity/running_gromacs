# Tiger (GPU)

The code can be built by following these commands:

```
$ ssh <YourNetID>@tiger3.princeton.edu
$ cd </path/to/your/software/directory>  # e.g., cd ~/software
$ wget https://raw.githubusercontent.com/PrincetonUniversity/running_gromacs/master/02_installation/ins/tigerGpu/tigerGpu.sh
# make modifications to tigerGpu.sh if needed (e.g., choose a different version)
$ bash tigerGpu.sh | tee build.log
```

The `tee` command is used to have the output go to the terminal as well as a file.
