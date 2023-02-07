# Build image
docker build -f Dockerfile --tag=ccbr_wgs_base:v0.1.0 .

# Tag image with version and reset latest
docker tag ccbr_wgs_base:v0.1.0 dnousome/ccbr_wgs_base:v0.1.0
docker tag ccbr_wgs_base:v0.1.0 dnousome/ccbr_wgs_base
#docker tag ccbr_wgs_base:v0.1.0 nciccbr/ccbr_wgs_base:v0.1.0
#docker tag ccbr_wgs_base:v0.1.0 nciccbr/ccbr_wgs_base

# Push image to DockerHub
docker push dnousome/ccbr_wgs_base:v0.1.0
docker push dnousome/ccbr_wgs_base:latest
#docker push nciccbr/ccbr_wgs_base:v0.1.0
#docker push nciccbr/ccbr_wgs_base:latest