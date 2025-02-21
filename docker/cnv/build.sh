##Build cnv
docker build --platform=linux/amd64 --tag ccbr_logan_cnv:v0.0.1 -f Dockerfile .

docker tag ccbr_logan_cnv:v0.0.1 dnousome/ccbr_logan_cnv:v0.0.1
docker tag ccbr_logan_cnv:v0.0.1 dnousome/ccbr_logan_cnv:latest

docker push dnousome/ccbr_logan_cnv:v0.0.1
docker push dnousome/ccbr_logan_cnv:latest


#singularity pull dnousome-ccbr_logan_cnv-v0.0.1.img docker://dnousome/ccbr_logan_cnv:v0.0.1
#docker run -it ccbr_logan_cnv:v0.0.1 
