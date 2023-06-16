for ENV_NAME in $(ls ../environments/ | grep ".yml" | grep -v "VEBA-preprocess_env"); 
do 
    ENV_NAME=$(echo $ENV_NAME | cut -f1 -d ".")
    echo $ENV_NAME
    
    time(bash build_docker_image.sh ${ENV_NAME})
done