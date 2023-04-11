# Docker

To build the image :

```shell
docker build -t <tag_name:version> -f <Dockerfile> .
```

To save the compressed image :

```shell
docker save <tag_name:version> | gzip > <filename.tar.gz>
```

To save it faster (with parallel gzip)

```shell
docker save <tag_name:version> | pigz --fast > <filename.tar.gz>
```