# VMD on the HPC Clusters

For intensive visualization work you can install [VMD](https://www.ks.uiuc.edu/Research/vmd/) in your `/home` directory.

## Installing VMD on Della

The following procedure can be used to install VMD on Della (and the other clusters):

```
$ ssh <YourNetID>@della-vis1.princeton.edu
$ cd software  # or another directory
$ wget https://www.ks.uiuc.edu/Research/vmd/vmd-1.9.4/files/alpha/vmd-1.9.4a57.bin.LINUXAMD64-CUDA102-OptiX650-OSPRay185.opengl.tar.gz
$ tar zxvf vmd-1.9.4a57.bin.LINUXAMD64-CUDA102-OptiX650-OSPRay185.opengl.tar.gz
$ cd vmd-1.9.4a57
```

Read the README file for the quick install directions then modify the `configure` file using a text editor (e.g., vim, emacs) with something like (replace <YourNetID> with your NetID):

```
$install_bin_dir="/home/<YourNetID>/software/vmd/bin";
$install_library_dir="/home/<YourNetID>/software/vmd/lib/$install_name";
```

Continue by running these commands:

```
$ ./configure
$ cd src
$ make install
# see "Updating your PATH" below
```

Use [MyDella](https://mydella.princeton.edu/) (or [MyStellar](https://mystellar.princeton.edu/) or [MyAdroit](https://myadroit.princeton.edu/)) to launch a graphical desktop by choosing "Interactive Apps" then "Desktop on Della Vis Nodes". When the session starts, click on the black terminal icon next to FireFox to launch a terminal. In the terminal, launch `vmd` with:

```
$ /home/<YourNetID>/software/vmd/bin/vmd
```

On tigressdata use:

```
$ vglrun /home/<YourNetID>/software/vmd/bin/vmd
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

### Working with Graphics

You should use [MyDella](https://mydella.princeton.edu/) or it is recommended that you install [TurboVNC](https://researchcomputing.princeton.edu/turbovnc) on your local machine (e.g., laptop) to work with VMD. If you decide to rely on X11 forwarding then make sure you are aware of [this page](https://researchcomputing.princeton.edu/sshX).

### More on TurboVNC

To use VMD on the head node of a cluster (please use MyDella/MyStellar/MyAdroit instead), first connect to `tigressdata.princeton.edu` using [TurboVNC](https://researchcomputing.princeton.edu/turbovnc) and then `ssh -X` to the desired cluster from tigressdata in a terminal. It will take you a few minutes to install and configure TurboVNC but you will find that the performance is excellent compared to X11 forwarding approaches.

Within TurboVNC, use these keyboard sequences for copy and paste:

```
Ctrl + Shift + C
Ctrl + Shift + V
```
