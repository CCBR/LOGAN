#Build.sh
docker build --platform=linux/amd64 --tag ccbr_logan_sv:v0.0.1 -f Dockerfile . 

docker tag ccbr_logan_sv:v0.0.1 dnousome/ccbr_logan_sv:v0.0.1
docker tag ccbr_logan_sv:v0.0.1 dnousome/ccbr_logan_sv:latest

docker push dnousome/ccbr_logan_sv:v0.0.1
docker push dnousome/ccbr_logan_sv:latest

##
#docker run -it ccbr_logan_sv:v0.0.1