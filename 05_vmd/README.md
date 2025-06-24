# VMD on the HPC Clusters

For intensive visualization work you can install [VMD](https://www.ks.uiuc.edu/Research/vmd/) in your `/home` directory.

## Installing VMD on Della

The following procedure can be used to install VMD on Della (and the other clusters):

```
$ ssh <YourNetID>@della-vis1.princeton.edu
$ cd software  # or another directory
$ wget https://www.ks.uiuc.edu/Research/vmd/alpha/vmd-2.0.0a5.bin.LINUXAMD64.tar.gz
$ tar zxvf vmd-2.0.0a5.bin.LINUXAMD64.tar.gz 
$ cd vmd-2.0.0a5
```

Read the README file for the quick install directions then modify the `configure` file using a text editor (e.g., vim, emacs) with something like (replace <YourNetID> with your NetID):

```
$install_bin_dir="/home/<YourNetID>/software/vmd/bin";
$install_library_dir="/home/<YourNetID>/software/vmd/lib/$install_name";
```

Continue by running these commands:

```
$ ./configure
using configure.options: LINUXAMD64 OPENGL OPENGLPBUFFER OPTIXRTRT FLTK TK ACTC CUDA CXX11 IMD LIBSBALL XINERAMA XINPUT LIBOPTIX LIBTACHYON LIBPNG ZLIB VRPN NETCDF COLVARS TCL PYTHON PTHREADS NUMPY SILENT
$ cd src
$ make install
Make sure /home/<YourNetID>/software/vmd/bin/vmd is in your path.
VMD installation complete.  Enjoy!

# see "Updating your PATH" below
```

Use [MyDella](https://mydella.princeton.edu/) (or [MyTiger](https://mytiger.princeton.edu) or [MyStellar](https://mystellar.princeton.edu/) or [MyAdroit](https://myadroit.princeton.edu/)) to launch a graphical desktop by choosing "Interactive Apps" then "Desktop on Della Vis Nodes". When the session starts, click on the black terminal icon next to FireFox to launch a terminal. In the terminal, launch `vmd` with:

```
$ /home/<YourNetID>/software/vmd/bin/vmd
```

In some cases you may benefit from running `vmd` under `vglrun`.

### Updating your PATH

Add VMD to your PATH environment variable:

```
$ vim ~/.bashrc  # or another text editor
```

Then add this line:

```
export PATH=$PATH:</path/to/vmd>  # e.g., export PATH=$PATH:/home/<YourNetID>/software/vmd/bin
```

After creating a new shell or sourcing your mofified `.bashrc` file, you can use `vmd`:

```
$ vmd <myfile.gro>
```
