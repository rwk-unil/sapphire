# DNANexus Applications

This directory holds DNANexus applications.

The applications can be built with :

```shell
dx build <path/to/app>
dx cp <name-of-applet> <Project:/path/where/to/store/the/applet>
```

Note that certain applications require binary resources and these must be copied beforehand.

There is one application that builds the other applications. `pp-toolkit-builder`, it removes the need to build them manually, it also takes care of copying the required binary files. Use this app to generate the other apps.