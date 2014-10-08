How to contribute to bdg-recipes
=========================

Thank you for sharing your code with the big data genomics project. We appreciate your contribution!

## Join the mailing list and our IRC channel

If you're not already on the ADAM developers list, [take a minute to join](http://bigdatagenomics.github.io/mail/).
It would be great if you'd introduce yourself to the group but it's not required. You can just
let your code do the talking for you if you like.

You can find us on Freenode IRC in the #adamdev room.

## Check the issue tracker

Before you write too much code, check the [open issues in the issue tracker](https://github.com/bigdatagenomics/bdg-recipes/issues?state=open)
to see if someone else has already filed an issue related to your work or is already working on it. If not, go ahead and 
[open a new issue](https://github.com/bigdatagenomics/bdg-recipes/issues/new).

## Submit your pull request

Github provides a nice [overview on how to create a pull request](https://help.github.com/articles/creating-a-pull-request).

Some general rules to follow:

* The recipes project is Python based, and uses the [Fabric](http://www.fabfile.org/) libraries and the [PEP-8](http://legacy.python.org/dev/peps/pep-0008/) styleguide.
* Do your work in [a fork](https://help.github.com/articles/fork-a-repo) of the repo.
* Create a branch for each feature/bug in that you're working on. These branches are often called "feature"
or "topic" branches.
* Use your feature branch in the pull request. Any changes that you push to your feature branch will automatically
be shown in the pull request.  If your feature branch is not based off the latest master, you will be asked to rebase
it before it is merged. This ensures that the commit history is linear, which makes the commit history easier to read.
* Before contributing code, check the [Github issue tracker](https://github.com/bigdatagenomics/bdg-recipes/issues).
If there is not an open ticket for what you would like to work on, please open it. When you submit your changes,
reference the issue from the pull request so that it will [be closed when your pull request is merged](https://github.com/blog/1506-closing-issues-via-pull-requests).
Also, please reference the issue number in your commit.
* Keep your pull requests as small as possible. Large pull requests are hard to review. Try to break up your changes
into self-contained and incremental pull requests, if need be, and reference dependent pull requests, e.g. "This pull
request builds on request #92. Please review #92 first."
* The first line of commit messages should start by referencing the issue number they fix (i.e., "[bdg-recipes-1]" indicates that
this commit fixes issue #1), followed by a short (<80 character) summary, followed by an empty line and then,
optionally, any details that you want to share about the commit.
