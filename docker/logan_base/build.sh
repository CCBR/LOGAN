
#Build image
docker build --platform linux/amd64 --tag ccbr_logan_base:v0.3.9 -f Dockerfile . 

docker tag ccbr_logan_base:v0.3.9 dnousome/ccbr_logan_base:v0.3.9
docker tag ccbr_logan_base:v0.3.9 dnousome/ccbr_logan_base:latest

docker push dnousome/ccbr_logan_base:v0.3.9
docker push dnousome/ccbr_logan_base:latest

#docker run -it ccbr_logan_base:v0.3.9

#Pull to CCBR
cd /data/CCBR_Pipeliner/SIFS
singularity pull dnousome-ccbr_logan_base-v0.3.9.img docker://dnousome/ccbr_logan_base:v0.3.9
