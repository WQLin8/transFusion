# transFusion: A Web-based Interactive Integrator for Spatial and Single-Cell Transcriptomes

<p align="center">
  <img src="Graphical_abstract.jpg" alt="transFusion Overview" width="150" align="right">
</p>

**transFusion** is a user-friendly, freely accessible tool developed to streamline the analysis of human and mouse scRNA-seq and ST data, leveraging advanced analytical techniques. It empowers researchers to gain deeper insights into spatially driven biological processes, including complex cell-cell communication networks, spatial regulatory networks, and tissue-specific functional organization.




## How to Use transFusion

You can use transFusion in **two ways**:

### 1. Run Locally with Docker
To run transFusion on your local computer, pull the image made in this repository from Docker Hub, use

```bash
sudo docker pull wlin8/transfusion:v1
```
the relevant Docker Hub repository can be found at <https://hub.docker.com/repository/docker/wlin8/transfusion>.

To use the platform, run

```bash
sudo docker run --rm --gpus all wlin8/transfusion:v1
```
Then, open **http://localhost:8080** in your browser.

### 2. Access via Cloud Server
Alternatively, you can directly access transFusion via the cloud:

➡️ **[Launch transFusion on Cloud](https://rhino-neat-woodcock.ngrok-free.app/app/transfusion)**  

## Citation

Weiqiang Lin, Xinyi Xiao, Hui Shen, Hongwen Deng*. transFusion: A Web-based Interactive Integrator for Spatial and Single-Cell Transcriptomes.

## Contact

For any questions and suggestions, please contact:

📧**Weiqiang Lin** (wlin8@tulane.edu)
