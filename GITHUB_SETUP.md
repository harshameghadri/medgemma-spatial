# GitHub Repository Setup Guide

**Status**: Ready to push ✅
**First Commit**: Created (966cf78)

---

## 📋 Pre-Push Checklist

- [x] Git repository initialized
- [x] .gitignore configured
- [x] Folder structure created
- [x] Initial commit made
- [ ] GitHub repository created
- [ ] Remote repository added
- [ ] Initial push completed
- [ ] Branch protection configured
- [ ] Secrets configured (for CI/CD)

---

## 🚀 Step-by-Step GitHub Setup

### Step 1: Create GitHub Repository

1. Go to https://github.com/new
2. Fill in repository details:
   ```
   Repository name: medgemma-spatial
   Description: Spatial transcriptomics analysis with MedGemma foundation models
   Visibility: Public (for portfolio) or Private

   ⚠️ DO NOT initialize with README (we already have one)
   ⚠️ DO NOT add .gitignore (we already have one)
   ⚠️ DO NOT add license (we already have MIT)
   ```
3. Click "Create repository"

### Step 2: Add Remote and Push

```bash
# Add GitHub as remote origin
git remote add origin https://github.com/YOUR_USERNAME/medgemma-spatial.git

# Verify remote
git remote -v

# Push initial commit
git push -u origin main

# Verify on GitHub
open https://github.com/YOUR_USERNAME/medgemma-spatial
```

### Step 3: Configure Branch Protection

1. Go to repository Settings → Branches
2. Click "Add rule" for `main` branch
3. Enable:
   - [x] Require pull request reviews before merging
   - [x] Require status checks to pass before merging
   - [x] Require branches to be up to date before merging
   - [x] Include administrators (recommended)

### Step 4: Set Up GitHub Secrets (for CI/CD)

Go to Settings → Secrets and variables → Actions → New repository secret

**Required Secrets:**

1. `HF_TOKEN` (for HuggingFace Spaces deployment)
   ```
   Get from: https://huggingface.co/settings/tokens
   Scope: write
   ```

2. `HF_USERNAME` (your HuggingFace username)
   ```
   Example: your-hf-username
   ```

**Optional Secrets:**

3. `CODECOV_TOKEN` (for code coverage)
   ```
   Get from: https://codecov.io/gh/YOUR_USERNAME/medgemma-spatial
   ```

### Step 5: Enable GitHub Actions

1. Go to repository → Actions tab
2. Click "I understand my workflows, go ahead and enable them"
3. Verify CI workflow runs on next push

---

## 🔄 Daily Git Workflow

### Starting New Feature

```bash
# Create feature branch
git checkout -b feature/scanpy-baseline

# Make changes
# ... edit files ...

# Stage and commit
git add notebooks/01_scanpy_baseline.ipynb
git commit -m "feat: Add scanpy baseline analysis notebook"

# Push to GitHub
git push origin feature/scanpy-baseline
```

### Creating Pull Request

1. Go to repository on GitHub
2. Click "Compare & pull request"
3. Fill in PR details:
   ```markdown
   ## Description
   Adds baseline spatial analysis using Scanpy

   ## Changes
   - Load and QC 10x Visium breast cancer dataset
   - Spatial leiden clustering
   - Export features to JSON

   ## Testing
   - [x] Ran notebook end-to-end
   - [x] Output JSON validates
   - [x] MPS acceleration working

   ## Related Issues
   Closes #1
   ```
4. Request review (if working with team)
5. Merge after CI passes

### Merging and Cleanup

```bash
# After PR is merged on GitHub
git checkout main
git pull origin main

# Delete local feature branch
git branch -d feature/scanpy-baseline

# Delete remote feature branch (if not auto-deleted)
git push origin --delete feature/scanpy-baseline
```

---

## 🏷️ Tagging Releases

### Create Version Tag

```bash
# Tag milestone (e.g., Week 1 completion)
git tag -a v0.1.0 -m "Week 1: Scanpy baseline complete"
git push origin v0.1.0

# Tag major version (e.g., deployment ready)
git tag -a v1.0.0 -m "Production deployment ready"
git push origin v1.0.0
```

### Semantic Versioning

- `v0.1.0`: Week 1 milestone
- `v0.2.0`: Week 2 milestone
- `v0.3.0`: Week 3 milestone (deployment)
- `v1.0.0`: Production ready (Week 4)

---

## 📊 GitHub Repository Settings

### About Section

```
Description: Spatial transcriptomics analysis with MedGemma foundation models
Website: https://huggingface.co/spaces/YOUR_USERNAME/medgemma-spatial
Topics: spatial-transcriptomics, scanpy, squidpy, medgemma, pathology,
        machine-learning, bioinformatics, streamlit
```

### Social Preview

Create a repository social image:
- Size: 1280x640px
- Content: Project logo + title + brief description

### Projects (Optional)

Create GitHub Project board for tracking:
- [ ] Week 1: Exploration
- [ ] Week 2: Integration
- [ ] Week 3: Deployment
- [ ] Week 4: Polish

---

## 🔐 Security Setup

### Enable Security Features

1. Settings → Security → Code security and analysis
2. Enable:
   - [x] Dependency graph
   - [x] Dependabot alerts
   - [x] Dependabot security updates
   - [x] Secret scanning

### Add SECURITY.md

```bash
# Create security policy
cat > SECURITY.md << 'EOF'
# Security Policy

## Supported Versions

| Version | Supported          |
| ------- | ------------------ |
| 1.0.x   | :white_check_mark: |
| < 1.0   | :x:                |

## Reporting a Vulnerability

Please report security vulnerabilities to: your.email@example.com
EOF

git add SECURITY.md
git commit -m "docs: Add security policy"
git push
```

---

## 🤖 Pre-commit Hook Setup

Install pre-commit hooks locally:

```bash
# Install pre-commit
pip install pre-commit

# Install git hooks
pre-commit install

# Test hooks
pre-commit run --all-files

# Now hooks run automatically on git commit
```

---

## 📝 Advanced Git Commands

### Squash Commits Before PR

```bash
# Squash last 3 commits
git rebase -i HEAD~3

# In editor, change "pick" to "squash" for commits to merge
# Save and exit

# Force push (only for feature branches, never main!)
git push --force-with-lease origin feature/your-branch
```

### Cherry-pick Commits

```bash
# Apply specific commit from another branch
git cherry-pick abc123
```

### Stash Changes

```bash
# Save work in progress
git stash

# Switch branches
git checkout other-branch

# Restore stashed changes
git stash pop
```

---

## 🎯 CI/CD Troubleshooting

### CI Fails on First Push

**Problem**: GitHub Actions might fail initially

**Solution**: Check these common issues:

1. **Secrets not set**: Add HF_TOKEN and HF_USERNAME
2. **Python version mismatch**: CI uses 3.10, 3.11
3. **Missing dependencies**: Update requirements.txt

### Pre-commit Hook Failures

**Problem**: Commits rejected by pre-commit

**Solution**:

```bash
# Fix formatting automatically
black src/ tests/
ruff check --fix src/ tests/

# Try commit again
git add -u
git commit -m "Your message"
```

---

## 🌐 HuggingFace Spaces Setup

### Create Space

1. Go to https://huggingface.co/new-space
2. Fill in:
   ```
   Name: medgemma-spatial
   License: MIT
   SDK: Streamlit
   Hardware: CPU Basic (free tier)
   Visibility: Public
   ```

### Link GitHub Repository

1. In Space settings → Files and versions
2. Link your GitHub repo: `YOUR_USERNAME/medgemma-spatial`
3. Enable automatic deployment on push to `main`

---

## 📦 Optional: Git LFS for Large Files

If you need to store files >100MB:

```bash
# Install Git LFS
brew install git-lfs  # macOS
# or: sudo apt-get install git-lfs  # Ubuntu

# Initialize
git lfs install

# Track large files
git lfs track "*.h5ad"
git lfs track "*.h5"

# Add .gitattributes
git add .gitattributes

# Commit as usual
git add data/large_file.h5ad
git commit -m "Add large dataset"
git push
```

**Warning**: GitHub has Git LFS bandwidth limits (1GB/month free)

---

## ✅ Verification Checklist

After setup, verify:

- [ ] Repository visible on GitHub
- [ ] CI workflow runs successfully
- [ ] Pre-commit hooks working locally
- [ ] Branch protection enabled
- [ ] Secrets configured
- [ ] README displays correctly
- [ ] License badge shows MIT
- [ ] Topics/tags added

---

## 🆘 Common Issues

### Push Rejected (Large Files)

```bash
# Error: file exceeds GitHub's 100MB limit

# Solution 1: Remove from history
git filter-branch --tree-filter 'rm -f path/to/large/file' HEAD
git push --force

# Solution 2: Use Git LFS (see above)
```

### Wrong Commit Message

```bash
# Amend last commit message (not pushed yet)
git commit --amend -m "New message"

# Already pushed? Create new commit
git commit --allow-empty -m "docs: Fix previous commit message"
```

### Accidentally Committed Secrets

```bash
# 1. Remove from latest commit
git rm --cached path/to/secret
git commit --amend -m "chore: Remove accidentally committed secret"

# 2. Remove from history
git filter-branch --force --index-filter \
  "git rm --cached --ignore-unmatch path/to/secret" \
  --prune-empty --tag-name-filter cat -- --all

# 3. Force push
git push --force

# 4. Rotate the exposed secret immediately!
```

---

## 📚 Resources

- [GitHub Docs](https://docs.github.com/)
- [Git Book](https://git-scm.com/book/en/v2)
- [GitHub Actions](https://docs.github.com/en/actions)
- [Pre-commit Hooks](https://pre-commit.com/)
- [Semantic Versioning](https://semver.org/)

---

## 🚀 Next Steps

After GitHub setup:

1. Create `develop` branch for active development
2. Set up GitHub Project board
3. Create first issue for Week 1 Day 1
4. Start development workflow

```bash
# Create develop branch
git checkout -b develop
git push -u origin develop

# Make develop the default branch on GitHub
# Settings → Branches → Default branch → develop
```

---

**Last Updated**: 2026-01-24
**Repository**: https://github.com/YOUR_USERNAME/medgemma-spatial (update after creation)
