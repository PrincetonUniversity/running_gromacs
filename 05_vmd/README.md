# VMD on the HPC Clusters

For intensive visualization work you can install [VMD](https://www.ks.uiuc.edu/Research/vmd/) in your `/home` directory on tigressdata. The `/scratch/gpfs` filesystems for each cluster and `/tigress` are
accessible from tigressdata. For example, for Perseus use the path: `/perseus/scratch/gpfs/<YourNetID>`.

For just doing quick checks of system configurations you can also install it in your `/home` directory on one of the
clusters (e.g., tiger).

## Working with Graphics

It is recommended that you install [TurboVNC](https://researchcomputing.princeton.edu/turbovnc) on your local machine (e.g., laptop) to work with VMD. If you decide to rely on X11 forwarding then make sure you are aware of [this page](https://researchcomputing.princeton.edu/sshX).

## Installing VMD on Tigressdata

The following procedure can be used to install VMD on tigressdata:

```
$ ssh -X <YourNetID>@tigressdata.princeton.edu
$ cd software  # or another directory
$ wget https://www.ks.uiuc.edu/Research/vmd/vmd-1.9.4/files/alpha/vmd-1.9.4a38.bin.LINUXAMD64-CUDA9-OptiX510-OSPRay170.opengl.tar.gz
$ tar zxvf vmd-1.9.4a38.bin.LINUXAMD64-CUDA9-OptiX510-OSPRay170.opengl.tar.gz
$ cd vmd-1.9.4a38
# read the README file for the quick install directions then
# modify the configure file with something like
$install_bin_dir="/home/<YourNetID>/software/vmd/bin";
$install_library_dir="/home/<YourNetID>/software/vmd/lib/$install_name";
$ cd src
$ make install
# launch VMD with the next command
$ vglrun /home/<YourNetID>/software/vmd/bin/vmd
# you coud add the above to your PATH in .bashrc
```

Note that `vglrun` is not necessary but on Tigressdata it allows for OpenGL rendering on the P100 GPU. On cluster head nodes it should be omitted.

### Updating your PATH

Add VMD to your PATH environment variable:

```
$ vim ~/.bashrc  # or another text editor
```

Then add this line:

```
export PATH=$PATH:</path/to/vmd>  # e.g., export PATH=$PATH:/home/<YourNetID>/software/vmd/bin
```

If you run VMD on a machine other than tigressdata then you should omit vglrun:

```
$ vmd <myfile.gro>
```

### More on TurboVNC

To use VMD on the head node of a cluster, first connect to `tigressdata.princeton.edu` using [TurboVNC](https://researchcomputing.princeton.edu/turbovnc) and then `ssh -X` to the desired cluster from tigressdata in a terminal. It will take you a few minutes to install and configure TurboVNC but you will find that the performance is excellent compared to X11 forwarding approaches.

Within TurboVNC, use these keyboard sequences for copy and paste:

```
Ctrl + Shift + C
Ctrl + Shift + V
```
