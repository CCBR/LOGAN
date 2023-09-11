# Build image
docker buildx build -f Dockerfile -t ccbr_wgs_base:v0.3.0 .

# Tag image with version and reset latest
docker tag ccbr_wgs_base:v0.3.0 dnousome/ccbr_wgs_base:v0.3.0
docker tag ccbr_wgs_base:v0.3.0 dnousome/ccbr_wgs_base
#docker tag ccbr_wgs_base:v0.1.0 nciccbr/ccbr_wgs_base:v0.1.0
#docker tag ccbr_wgs_base:v0.1.0 nciccbr/ccbr_wgs_base

# Push image to DockerHub
docker push dnousome/ccbr_wgs_base:v0.3.0 ##Push to v0.3.0 branch
docker push dnousome/ccbr_wgs_base:latest ##Push to Latest tag
#docker push nciccbr/ccbr_wgs_base:v0.1.0
#docker push nciccbr/ccbr_wgs_base:latest