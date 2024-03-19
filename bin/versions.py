#!/usr/bin/env python

# Stolen and modified from 
# https://github.com/nf-core/rnaseq/blob/b89fac32650aacc86fcda9ee77e00612a1d77066/modules/nf-core/custom/dumpsoftwareversions/templates/dumpsoftwareversions.py#L4

""" Reformat versions yml file for multiqc """

import yaml
from textwrap import dedent

def _make_versions_html(versions):
    """ Reformat versions yml file for multiqc """

    html = [
        dedent(
            """\\
            <style>
            #nf-core-versions tbody:nth-child(even) {
                background-color: #f2f2f2;
            }
            </style>
            <table class="table" style="width:100%" id="nf-core-versions">
                <thead>
                    <tr>
                        <th> Process Name </th>
                        <th> Software </th>
                        <th> Version  </th>
                    </tr>
                </thead>
            """
        )
    ]
    for process, tmp_versions in sorted(versions.items()):
        html.append("<tbody>")
        for i, (tool, version) in enumerate(sorted(tmp_versions.items())):
            html.append(
                dedent(
                    f"""\\
                    <tr>
                        <td><samp>{process if (i == 0) else ''}</samp></td>
                        <td><samp>{tool}</samp></td>
                        <td><samp>{version}</samp></td>
                    </tr>
                    """
                )
            )
        html.append("</tbody>")
    html.append("</table>")
    return "\n".join(html)


def main():
    """Load all version files and generate merged output."""

    with open("versions.yml") as f:
        versions_by_process = yaml.load(f, Loader=yaml.BaseLoader) 

    versions_by_module = {}
    for process, process_versions in versions_by_process.items():
        module = process.split(":")[-1]
        try:
            if versions_by_module[module] != process_versions:
                raise AssertionError(
                    "There's something wrong with the designated containers of this workflow"
                )
        except KeyError:
            versions_by_module[module] = process_versions

    versions_mqc = {
        "id": "software_versions",
        "section_name": "CECRET Software Versions",
        "section_href": "https://github.com/UPHL-BioNGS/Grandeur",
        "plot_type": "html",
        "description": "Collected at run time from the software output.",
        "data": _make_versions_html(versions_by_module),
    }

    with open("software_versions.yml", "w") as f:
        yaml.dump(versions_by_module, f, default_flow_style=False)
    with open("software_versions_mqc.yml", "w") as f:
        yaml.dump(versions_mqc, f, default_flow_style=False)

if __name__ == "__main__":
    main()