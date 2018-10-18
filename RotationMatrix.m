function [RotationX,RotationY,RotationZ] = RotationMatrix(EllipsoidAngleX,EllipsoidAngleY,EllipsoidAngleZ)



    RotationX = [   1             0                    0                0
                    0   cos(EllipsoidAngleX)  -sin(EllipsoidAngleX)     0
                    0   sin(EllipsoidAngleX)  cos(EllipsoidAngleX)      0
                    0             0                    0                1  ];
    
    RotationY = [   cos(EllipsoidAngleY)   0   sin(EllipsoidAngleY)     0
                             0             1           0                0
                    -sin(EllipsoidAngleY)  0  cos(EllipsoidAngleY)      0
                             0             0           0                1  ];
                         
    RotationZ = [   cos(EllipsoidAngleZ)   -sin(EllipsoidAngleZ)   0    0
                    sin(EllipsoidAngleZ)   cos(EllipsoidAngleZ)    0    0
                             0                      0              1    0
                             0                      0              0    1  ];