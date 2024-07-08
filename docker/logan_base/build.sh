
# Build image
#docker buildx create --platform linux/amd64 --use
#docker buildx use upbeat_ganguly
#docker buildx inspect upbeat_ganguly
#docker buildx build --platform linux/amd64 -f Dockerfile -t dnousome/ccbr_logan_base:v0.3.0 -t dnousome/ccbr_logan_base:latest --push .

docker build --platform linux/amd64 --tag ccbr_logan_base:v0.3.6 -f Dockerfile . 

docker tag ccbr_logan_base:v0.3.6 dnousome/ccbr_logan_base:v0.3.6
docker tag ccbr_logan_base:v0.3.6 dnousome/ccbr_logan_base


docker push dnousome/ccbr_logan_base:v0.3.6
docker push dnousome/ccbr_logan_base:latest




# Tag image with version and reset latest
#docker tag ccbr_wgs_base:v0.1.0 nciccbr/ccbr_wgs_base:v0.1.0
#docker tag ccbr_wgs_base:v0.1.0 nciccbr/ccbr_wgs_base

# Push image to DockerHub
#docker push nciccbr/ccbr_wgs_base:v0.1.0
#docker push nciccbr/ccbr_wgs_base:latest
