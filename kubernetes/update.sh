#!/usr/bin/env zsh

set -e
TAG=`cat ../VERSION`
echo $TAG
echo "Updating deployment..."
if [[ `uname` == 'Darwin' ]]
then
    # Mac uses BSD version of sed, which differs from GNU implementation. For
    # details see:
    # https://unix.stackexchange.com/questions/401905/bsd-sed-vs-gnu-sed-and-i
    # Idea for future improvement: parse `sed --version` for "GNU" or "BSD"
    # instead of evaluating `uname`.
    sed -i '' -e "s/\(image:.*:\).*$/\1$TAG/g" deployment.yaml
else
    sed -i -e "s/\(image:.*:\).*$/\1$TAG/g" deployment.yaml
fi
echo "Applying deployment..."
kubectl apply -f deployment.yaml
echo "Commiting..."
git add deployment.yaml
git commit -m "deployment $TAG"
echo "Done."
