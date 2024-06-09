for DOCKER_IMAGE in $(docker images --format "{{.Repository}}:{{.Tag}}" | grep "veba")
do
    echo $DOCKER_IMAGE
    docker push $DOCKER_IMAGE
    sleep 1
done
