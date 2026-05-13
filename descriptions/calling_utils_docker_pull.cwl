#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool

hints:
  - class: DockerRequirement
    dockerPull: ACCOUNT/calling_utils:VERSION

baseCommand: [echo, "pulled calling_utils docker image!"]

inputs: []

outputs:
  - id: output_log_txt
    type: stdout

stdout: log.txt

doc: |
    Empty CWL to be used for pre-pulling the calling_utils Docker image
