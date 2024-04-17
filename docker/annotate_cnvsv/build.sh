##BUILD cnv/sv
docker build --platform linux/amd64 --tag ccbr_annotate_cnvsv:v0.0.1 -f Dockerfile . 

docker tag ccbr_annotate_cnvsv:v0.0.1 dnousome/ccbr_annotate_cnvsv:v0.0.1
docker tag ccbr_annotate_cnvsv:v0.0.1 dnousome/ccbr_annotate_cnvsv

docker push dnousome/ccbr_annotate_cnvsv:v0.0.1
docker push dnousome/ccbr_annotate_cnvsv:latest


