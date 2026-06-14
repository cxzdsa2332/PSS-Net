---
name: commit
description: Auto-commit the current PSS-Net changes — inspect the working tree, group related changes, and create a git commit whose message summarizes the main content of this round of edits (Chinese summary + Co-Authored-By trailer). Use when the user asks to commit, 提交, or save changes to git.
---

# commit — PSS-Net 自动提交本次修改

检查当前工作区改动，归纳**本次修改的主要内容**，生成一条规范的 git commit。

## 执行步骤

1. **查看现状**（并行运行，单条消息内）：
   - `git status` —— 看暂存/未暂存/未跟踪文件。
   - `git diff` 与 `git diff --staged` —— 看实际改动内容。
   - `git log --oneline -5` —— 对齐已有提交信息的语言与风格。

2. **判断范围**：
   - 默认提交本次会话产生的、彼此相关的改动。若工作区还混有与本次任务无关的旧改动，先向用户确认是否一并提交，不要默认 `git add -A` 吞掉一切。
   - 不提交构建产物或被 `.gitignore` 忽略的文件（如 `manuscript/main.pdf`、LaTeX 中间件、`.DS_Store`）。

3. **暂存**：`git add` 明确列出要提交的文件（优于 `-A`）。

4. **写 commit message**（中文，与仓库历史一致）：
   - 首行：`<动词>: <一句话概括本次主要改动>`（≤ 50 字，如「新增 manuscript LaTeX 骨架与 commit skill」）。
   - 空行后用 `-` 列出 2–5 条主要变更点（按文件/模块归纳，说清「做了什么」而非逐行罗列）。
   - 结尾固定加 trailer：
     ```
     Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>
     ```
   - 用 HEREDOC 传入以保证多行格式：
     ```bash
     git commit -m "$(cat <<'EOF'
     <首行概括>

     - <要点1>
     - <要点2>

     Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>
     EOF
     )"
     ```

5. **校验**：`git status` 确认提交成功、工作区干净（或仅剩有意未提交的内容）；向用户回报提交哈希与一句话摘要。

## 约束
- **仅在用户要求时提交**；不主动 `git push`（除非用户明确要求）。
- 若当前在主分支且改动较大或属功能性变更，先提示是否应另开分支再提交。
- 不修改 git 历史（不 `commit --amend` / `rebase`，除非用户明确要求）。
- 提交信息只描述事实，不夸大；测试未跑就不要在信息里声称已验证。
