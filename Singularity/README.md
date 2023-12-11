# Singularity

## 1. Automatic creation of Singularity recipe

```shell
spython recipe ../Docker/Dockerfile > ./Dockerfile.srecipe
```

## 2. Manual Editing of the recipe

 - Exclude `rm -rf /tmp*`-like command for safety
 - Add the `STATIC_BINS=y` option for `make`
 - Add compile commands for `bin_tools` and `pp_update`
