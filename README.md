Check file changes (staged/unstaged) and number of commits ahead/behind upstream
```git
git status
```

Download changes on all branches and apply upstream changes to this branch
```git
git pull
```

Adding files to be committed (and variations) (-A for all changes)
```git
git add <filename>
git add -A
```

Committing changes
```git
git commit         // this launches vim where you can leave a message
git commit -m "message goes here"
```

Upload your commits... (you might have to set origin if it's the first commit on a branch - git gives you a guide)
```git
git push
```

Change to a different branch 'feature/something_existing'
```git
git checkout feature/something_existing
```
  
Create and change to new branch 'feature/something_new' and take current changes with you
```git
git checkout -b feature/something_new
```
  
merge changes from another branch i.e. 'master', 'feature/something_new' etc. into current branch
```git
git merge master
```
