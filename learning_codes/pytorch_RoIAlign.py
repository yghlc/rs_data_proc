#!/usr/bin/env python
# Filename: pytorch_RoIAlign.py 
"""
introduction:

authors: Huang Lingcao
email:huanglingcao@gmail.com
add time: 18 July, 2025
"""
import os.path

import torch
from clip import clip

from PIL import Image

data_dir = os.path.expanduser('~/codes/PycharmProjects/BigImageMapper/img_classification')

def compare_clip_model_visual_vs_encode_image():
    # Load CLIP model (ViT-B/32 backbone)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    model, preprocess = clip.load("ViT-B/32", device=device)
    # model, preprocess = clip.load("RN50", device=device)

    # Load and preprocess the image
    image_path = os.path.join(data_dir,"intersection88.tif")
    image = Image.open(image_path).convert("RGB")  # Ensure image is in RGB format
    image_tensor = preprocess(image).unsqueeze(0).to(device)  # Shape: (1, 3, H, W)

    # Extract Intermediate Visual Features using model.visual
    with torch.no_grad():  # Disable gradient computation
        visual_features = model.visual(image_tensor)

    print("Visual Features Shape (from model.visual):", visual_features.shape)
    # print(visual_features)
    # For ViT-B/32:
    # - Shape = (1, num_tokens, 512), where num_tokens = (H/32)*(W/32) + 1

    # Extract Final Image Embedding using model.encode_image
    with torch.no_grad():
        image_embedding = model.encode_image(image_tensor)

    print("Final Image Embedding Shape (from model.encode_image):", image_embedding.shape)
    # print(image_embedding)
    # For ViT-B/32:
    # - Shape = (1, 512), a single vector embedding for the entire image

def main():
    # device = "cuda" if torch.cuda.is_available() else "cpu"
    compare_clip_model_visual_vs_encode_image()


    pass


if __name__ == '__main__':
    main()
