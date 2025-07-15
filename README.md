# PoTMoP: PTM-Processing Pipeline 

## Run in docker / enroot 
This pipeline is part of a bigger pipeline that can be run with enroot on an HPC. For information about that, refer to: <TODO insert github link>. 

## Run localy
Create a virtual environment and install the requirements.txt. 
I'd recommend to create a `launch.json` in vscode: 

```
{
    "version": "0.2.0",
    "configurations": [       
        {
            "name": "Init Pipeline",
            "type": "python",
            "request": "launch",
            "program": "main.py",
            "args": ["init"],
            "console": "integratedTerminal"
        },
        {
            "name": "Run Pipeline - Prepare",
            "type": "python",
            "request": "launch",
            "program": "main.py",
            "args": ["prepare"],
            "console": "integratedTerminal"
        },
        {
            "name": "Run Pipeline - Create DB",
            "type": "python",
            "request": "launch",
            "program": "main.py",
            "args": ["create_db"],
            "console": "integratedTerminal"
        }, 
        {
            "name": "Run Pipeline - Filter sig. Proteins",
            "type": "python",
            "request": "launch",
            "program": "main.py",
            "args": ["filter_significant_proteins"],
            "console": "integratedTerminal"
        }, 
        {
            "name": "Run Pipeline - Process new Af2 Data",
            "type": "python",
            "request": "launch",
            "program": "main.py",
            "args": ["process_new_alphafold2_runs"],
            "console": "integratedTerminal"
        },
        {
            "name": "Run Pipeline - Extract and format PTMs",
            "type": "python",
            "request": "launch",
            "program": "main.py",
            "args": ["preprocess"],
            "console": "integratedTerminal"
        },
        {
            "name": "Run Pipeline - StructureMap preprocessing",
            "type": "python",
            "request": "launch",
            "program": "main.py",
            "args": ["structure_map_preprocessing"],
            "console": "integratedTerminal"
        },
        {
            "name": "Analysis",
            "type": "python",
            "request": "launch",
            "program": "first_analysis_python.py",
            "args": ["init"],
            "console": "integratedTerminal",
            "justMyCode": false,
        },
        {
            "name": "Python: pytest",
            "type": "python",
            "request": "launch",
            "program": "${workspaceFolder}/venv/bin/pytest", 
            "args": [
              "${workspaceFolder}/tests", 
              "--maxfail=1",
              "--disable-warnings"
            ],
            "env": {
              "PYTHONPATH": "${workspaceFolder}/src"
            },
            "justMyCode": false
          }
    ]
}
```