# Build image

docker build --platform linux/amd64 --tag ccbr_logan_base:v0.3.4 -f Dockerfile . 

docker tag ccbr_lofreq:v0.0.1 dnousome/ccbr_lofreq:v0.0.1
docker push dnousome/ccbr_lofreq:v0.0.1

docker push dnousome/ccbr_logan_base:latest



