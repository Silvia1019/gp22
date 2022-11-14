# Assignment 5

Name: Haitong Shi

Legi-Nr: 20-960-340

## Required results
Edit this 'README.md' file to report all your results. You only need to update the tables in the reports section by adding screenshots and reporting results.

### Tasks

1) **Multiresolution Mesh Editing**: Provide screenshots for 4 different deformed meshes (woody-lo, woody-hi, hand, and cylinder). For each example, provide a rendering of S, B, B' and S'. (questions 1.1 - 1.4 of the assignment sheet)

2) **Real time mesh editing**: Provide animated gifs or short videos for 4 different deformed meshes (bar, bumpy_plane, camel_head, and cactus) showing that your algorithm can perform in real time. (question 1.5 of the assignment sheet)

3) **Deformation transfer**: Discuss and show the differences to the results obtained with the high-frequency detail transfer from part 1.4. on 4 different meshes (bar, bumpy_plane, camel_head, and cactus). (part 2 of the assignment sheet)



## Reports
### 1.1 - 1.4 Multiresolution Mesh Editing
| model name     | S     |  B    |  B'   |  S'   |
| :-----------:  | ----- | ----- | ----- | ----- |
| woody-lo       |<img align="center" src="./res/woody_lo_s.png" width="300">| <img align="center"  src="./res/woody_lo_b.png" width="300"> |<img align="center" src="./res/woody_lo_bp.png" width="300">| <img align="center"  src="./res/woody_lo_sp.png" width="300"> |
| woody-hi       |<img align="center" src="./res/woody_hi_s.png" width="300">| <img align="center"  src="./res/woody_hi_b.png" width="300"> |<img align="center" src="./res/woody_hi_bp.png" width="300">| <img align="center"  src="./res/woody_hi_sp.png" width="300"> |
| hand           |<img align="center" src="./res/hand_s.png" width="300">| <img align="center"  src="./res/hand_b.png" width="300"> |<img align="center" src="./res/hand_bp.png" width="300">| <img align="center"  src="./res/hand_sp.png" width="300"> |
| cylinder       |<img align="center" src="./res/cylin_s.png" width="300">| <img align="center"  src="./res/cylin_b.png" width="300"> |<img align="center" src="./res/cylin_bp.png" width="300">| <img align="center"  src="./res/cylin_sp.png" width="300"> |

### 1.5 Real time mesh editing

Show real time mesh editing using animated gifs or short videos. *Max 15 seconds per gif, better if 5 to 10 seconds*.

(The original video is less than 15 seconds but gifs are a bit slower and longer. Here are links to videos.)

[bar](./res/bar.mov) 

[bumpy](./res/bumpy.mov) 

[camel](./res/camel.mov) 

[cactus](./res/cac.mov)

| model name     |   S' - real time   |
| :-----------:  | -----  |
| bar            | <img align="center"  src="./res/bar.gif" width="300"> |
| bumpy_plane    | <img align="center"  src="./res/bumpy.gif" width="300"> |
| camel_head     | <img align="center"  src="./res/camel.gif" width="300"> |
| cactus         | <img align="center"  src="./res/cac.gif" width="300"> |


### 2 Deformation transfer
| model name     | High-freq detail transfer             | Deformation transfer                 |
| :-----------:  | ------------------------------------- |------------------------------------- |
| bar            |<img align="center" src="./res/bar_multi.png" width="300">| <img align="center"  src="./res/bar_de.png" width="300"> |
| bumpy_plane    |<img align="center" src="./res/bumpy_multi.png" width="300">| <img align="center"  src="./res/bumpy_de.png" width="300"> |
| camel_head     |<img align="center" src="./res/camel_multi.png" width="300">| <img align="center"  src="./res/camel_de.png" width="300"> |
| cactus         |<img align="center" src="./res/cac_multi.png" width="300">| <img align="center"  src="./res/cac_de.png" width="300"> |


#### Observations

|      | High-freq detail transfer             | Deformation transfer                 |
| :-----------:  | ------------------------------------- |------------------------------------- |
| Your Comments  | One limitation of this approach is that the difference between S and B must be sufficiently small, such that S can be represented as a height field over B. The resulting displaced surface may have local self-intersections when the curvature of the deformed base surface B′ is too high in relation to the displacement length. We can see in bumpy_plane, local details are changed and the shape is distorted.  | This method combines the advantages of differential coordinates and explicit multiresolution decomposition. The rotation is smoothier and the shape and details are well preserved especially in bumpy_plane. But due to larger matrix multiplication, this method is slower than the former in my implementation. |
 