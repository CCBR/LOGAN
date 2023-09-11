# Build image
docker buildx create --platform linux/amd64 --use
docker buildx build -f Dockerfile -t dnousome/ccbr_logan_base:v0.3.0 -t dnousome/ccbr_logan_base:latest  --platform=linux/amd64 --push .

# Tag image with version and reset latest
#docker tag ccbr_wgs_base:v0.1.0 nciccbr/ccbr_wgs_base:v0.1.0
#docker tag ccbr_wgs_base:v0.1.0 nciccbr/ccbr_wgs_base

# Push image to DockerHub
#docker push nciccbr/ccbr_wgs_base:v0.1.0
#docker push nciccbr/ccbr_wgs_base:latest