<!-- badges: start -->
[![R-CMD-check](https://github.com/yerkes-gencore/gencoreSC/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yerkes-gencore/gencoreSC/actions/workflows/R-CMD-check.yaml)
[![Docker Pulls](https://img.shields.io/docker/pulls/yerkesgencore/gencore-singlecell-rstudio)](https://hub.docker.com/r/yerkesgencore/gencore-singlecell-rstudio)
<!-- badges: end -->

# gencoreSC - Utility functions for single cell analyses

## This is an R package

This is intended as a living codebase for useful custom R functions for our single cell analysis workflows downstream of generating UMI count tables (e.g. via cellranger). 

It's structured as an R package to keep things tidy, well documented and version controlled. 

# Installation

You can install this package using `devtools` or `remotes`.

```
devtools::install_github('yerkes-gencore/gencoreSC')
```

Some dependencies don't automatically install for various reasons. If you have issues installing the package due to missing dependencies, try installing them directly,
then reinstalling the package. For example, the `S4Vectors` seems to throw issues occasionally. Run `BiocManager::install("S4Vectors")` then reinstall gencoreSC. 

The package has been tested to install with the docker image

https://hub.docker.com/r/yerkesgencore/gencore-singlecell-rstudio

# Making changes

Before making changes, please also review tutorials on developing/maintaining simple R packages, such as this one: https://kbroman.org/pkg_primer/. The essential stuff is maybe a 30 min read and if we all familiarize with the basics, then it shouldn't be too hard to keep this simple and useful.

If you want to make any changes, follow this workflow

1. clone this repo into an isolated working directory (we suggest `.../illumina/runs/analyst/<your_name>`)
2. Fetch updates from the github repo via `git pull`
3. Create a new branch from the `devel` branch via `git branch <new_branch> origin/devel`
4. Make changes on the new branch. When you're finished making changes, be sure to run 
`devtools::check()` to ensure the documentation is updated and their are no
major issues. When you're satisfied with the changes, push the changes to Github.
5. Open a pull request to merge your new branch into the devel branch. Ideally someone else should review the changes before merging, and the RMD Check action should pass.
6. Once the feature has been merged into the devel branch, you should safely delete that branch with `git branch -d <your-branch>`. You can always recover it with `git checkout <your-branch> <sha>`, where `<sha>` is the identifying SHA string for the commit at the tip of that branch (you can always find that in your git history).
7. Once sufficient changes have been made to the devel branch to prompt an update to the package, modify the package description to update the version
8. Release a new version of the package

See this writeup on the git flow workflow 

https://medium.com/@patrickporto/4-branching-workflows-for-git-30d0aaee7bf

We could switch to other workflows as well, but for now we'll use git flow since we don't really need continuous deployment

https://rakeshjain-devops.medium.com/fix-to-tip-of-your-current-branch-is-behind-its-remote-counterpart-git-error-eb75f719c2d5
