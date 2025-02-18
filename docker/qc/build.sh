#Build.sh
docker build --platform=linux/amd64 --tag ccbr_logan_qc:v0.0.2 -f Dockerfile . 

docker tag ccbr_logan_qc:v0.0.2 dnousome/ccbr_logan_qc:v0.0.2
docker tag ccbr_logan_qc:v0.0.2 dnousome/ccbr_logan_qc:latest

docker push dnousome/ccbr_logan_qc:v0.0.2
docker push dnousome/ccbr_logan_qc:latest

##
#docker run -it ccbr_logan_qc:v0.0.2