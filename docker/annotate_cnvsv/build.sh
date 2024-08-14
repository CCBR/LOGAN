##BUILD cnv/sv
docker build --platform linux/amd64 --tag ccbr_annotate_cnvsv:v0.0.2 -f Dockerfile . 

docker tag ccbr_annotate_cnvsv:v0.0.2 dnousome/ccbr_annotate_cnvsv:v0.0.2
docker tag ccbr_annotate_cnvsv:v0.0.2 dnousome/ccbr_annotate_cnvsv:latest

docker push dnousome/ccbr_annotate_cnvsv:v0.0.2
docker push dnousome/ccbr_annotate_cnvsv:latest


