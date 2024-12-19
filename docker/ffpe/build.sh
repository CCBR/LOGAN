#Build.sh
docker build --platform=linux/amd64 --tag ccbr_logan_ffpe:v0.0.1 -f Dockerfile . 

#Test Docker Build
#docker run -it ccbr_logan_ffpe:v0.0.1
#

docker tag ccbr_logan_ffpe:v0.0.1 dnousome/ccbr_logan_ffpe:v0.0.1
docker tag ccbr_logan_ffpe:v0.0.1 dnousome/ccbr_logan_ffpe:latest

docker push dnousome/ccbr_logan_ffpe:v0.0.1
docker push dnousome/ccbr_logan_ffpe:latest

