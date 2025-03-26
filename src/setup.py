from setuptools import setup, find_packages

setup(
    name="plmpg",
    version="0.1.0",
    description="Unified ESM2 + ESM3 inference and processing tools",
    author="Guy 'Wayyne' Dayhoff",
    author_email="gdayhoff@rx.umaryland.edu",
    packages=find_packages(where=".", include=["plmpg", "plmpg.*"]),
    package_dir={"": "."},
    python_requires=">=3.10",
    install_requires=[],  # Already installed from KaML-ESM_env.txt
    entry_points={
        "console_scripts": [
            # Uncomment to expose main.py as a CLI
            # "kaml-main=plmpg.main:main",
        ]
    },
    include_package_data=True,
    zip_safe=False,
)
