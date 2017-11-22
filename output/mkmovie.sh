#!/bin/bash
ffmpeg -r 6 -i image_%04d.bmp -vcodec libx264 -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -pix_fmt yuv420p -r 60 out.mp4 &> output_log.txt
