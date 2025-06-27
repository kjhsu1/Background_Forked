#!/usr/bin/env bash
# Script will complement the codex "Setup Script"
# In AGENTS.md, tell it to source this before anything

INSTALL_DIR="$HOME/.micromamba"
BIN_DIR="$INSTALL_DIR/bin"

export PATH="$BIN_DIR:$PATH"

micromamba activate codex-env