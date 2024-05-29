# Build image

docker build --platform linux/amd64 --tag ccbr_logan_hmftools:v0.0.1 -f Dockerfile . 

docker tag ccbr_logan_hmftools:v0.0.1 dnousome/ccbr_logan_hmftools:v0.0.1
docker push dnousome/ccbr_logan_hmftools:v0.0.1

docker push dnousome/ccbr_logan_hmftools:latest



