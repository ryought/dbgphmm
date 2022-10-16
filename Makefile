build:
	maturin build && pip install --force-reinstall target/wheels/dbgphmm-0.1.0-cp310-cp310-macosx_11_0_arm64.whl
