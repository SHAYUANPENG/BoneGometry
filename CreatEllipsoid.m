function [P,PXYZ,centroid,axisLength,angles,F,V,label] = CreatEllipsoid

    global Ellipsoids
    
    Ellipsoid_precision = Ellipsoids.Ellipsoid_precision ;
    
    EllipsoidAxisXmin =  Ellipsoids.EllipsoidAxisXmin;
    EllipsoidAxisXmax =  Ellipsoids.EllipsoidAxisXmax;
    EllipsoidAxisYmin =  Ellipsoids.EllipsoidAxisYmin;
    EllipsoidAxisYmax =  Ellipsoids.EllipsoidAxisYmax;
    EllipsoidAxisZmin =  Ellipsoids.EllipsoidAxisZmin;
    EllipsoidAxisZmax =  Ellipsoids.EllipsoidAxisZmax;    %椭球X、Y、Z三轴最小最大长度，确定随机变量范围
    
    EllipsoidAngleXmin =  Ellipsoids.EllipsoidAngleXmin;
    EllipsoidAngleXmax =  Ellipsoids.EllipsoidAngleXmax;
    EllipsoidAngleYmin =  Ellipsoids.EllipsoidAngleYmin;
    EllipsoidAngleYmax =  Ellipsoids.EllipsoidAngleYmax;
%    EllipsoidAngleZmin =  Ellipsoids.EllipsoidAngleZmin;
%    EllipsoidAngleZmax =  Ellipsoids.EllipsoidAngleZmax;          %椭球与X、Y、Z三轴夹角最小最大长度，确定随机变量范围

    EllipsoidDistanceXmin = Ellipsoids.EllipsoidDistanceXmin;
    EllipsoidDistanceXmax = Ellipsoids.EllipsoidDistanceXmax;
    EllipsoidDistanceYmin = Ellipsoids.EllipsoidDistanceYmin;
    EllipsoidDistanceYmax = Ellipsoids.EllipsoidDistanceYmax;
    EllipsoidDistanceZmin = Ellipsoids.EllipsoidDistanceZmin;
    EllipsoidDistanceZmax = Ellipsoids.EllipsoidDistanceZmax;       %两个椭球之间的最小距离
    
    Height_min = Ellipsoids.Height_min;
    Height_max = Ellipsoids.Height_max;   %椭球中心点出现的区域
    
    Thickness1 = Ellipsoids.Thickness1;
    Thickness2 = Ellipsoids.Thickness2;
    
    p = Ellipsoids.p;
    p_angle = Ellipsoids.p_angle.matrix;
    
     % 根据0.7和0.3的概率，使椭球落在层内或者层间
    a = rand(1,1);
    if a>p
        label = 1; %label为1代表层间椭球
    else
        label = 0; %label为0代表层内椭球
    end
    % 控制落在各圆环的中心线上
    innerRadius = Ellipsoids.BoneInnerRadius;
    if label == 0
      %假设最里面一层和最外面一层层内没有椭球，共有八条圆环中心线，若有椭球，可以添加两条中心线
      Ellipsoidradius1 = [innerRadius +Thickness1+Thickness2/2, innerRadius+Thickness2+3/2*Thickness1,innerRadius+2*Thickness1+3/2*Thickness2,...
                         innerRadius+5/2*Thickness1+2*Thickness2,innerRadius+3*Thickness1+5/2*Thickness2,innerRadius+7/2*Thickness1+3*Thickness2,...
                         innerRadius+4*Thickness1+7/2*Thickness2,innerRadius+9/2*Thickness1+4*Thickness2];
      C=randperm(numel(Ellipsoidradius1));
      Ellipsoidradius=Ellipsoidradius1(C(1));       %随机抽取一个半径
    else
        Ellipsoidradius1 = [innerRadius +Thickness1,innerRadius +Thickness1+Thickness2,innerRadius +2*Thickness1+Thickness2,...
                            innerRadius +2*Thickness1+2*Thickness2,innerRadius +3*Thickness1+2*Thickness2,innerRadius +3*Thickness1+3*Thickness2,...
                            innerRadius +4*Thickness1+3*Thickness2,innerRadius +4*Thickness1+4*Thickness2,innerRadius +5*Thickness1+4*Thickness2];
        C=randperm(numel(Ellipsoidradius1));
        Ellipsoidradius=Ellipsoidradius1(C(1));       %随机抽取一个半径
    end
    EllipsoidHeight = Height_min + (Height_max-Height_min).*rand(1,1);
    EllipsoidTheta = 360*rand(1,1)/180*pi;
    
    EllipsoidX = Ellipsoidradius*cos(EllipsoidTheta);
    EllipsoidY = Ellipsoidradius*sin(EllipsoidTheta);
    EllipsoidZ = EllipsoidHeight;
    centroid = [EllipsoidX,EllipsoidY,EllipsoidZ];
    
    %随机椭球三个轴长度
    EllipsoidAxisX = EllipsoidAxisXmin + (EllipsoidAxisXmax-EllipsoidAxisXmin).*rand(1,1);
    EllipsoidAxisY = EllipsoidAxisYmin + (EllipsoidAxisYmax-EllipsoidAxisYmin).*rand(1,1);
    EllipsoidAxisZ = EllipsoidAxisZmin + (EllipsoidAxisZmax-EllipsoidAxisZmin).*rand(1,1);
    axisLength = [EllipsoidAxisX,EllipsoidAxisY,EllipsoidAxisZ];
    
    %随机椭球三个角度,根据概率控制是否发生旋转
    label_angle = rand(1,1);
    if  label_angle>0  && label_angle<p_angle(1)
        EllipsoidAngleX = 0;
        EllipsoidAngleY = 0;
        EllipsoidAngleZ = 0;
    elseif  label_angle>p_angle(1) && label_angle<(p_angle(1)+p_angle(2))
        EllipsoidAngleX = EllipsoidAngleXmin + (EllipsoidAngleXmax-EllipsoidAngleXmin).*rand(1,1);
        EllipsoidAngleY = 0;
        EllipsoidAngleZ = 0;
    elseif  label_angle>(p_angle(1)+p_angle(2)) && label_angle<(p_angle(1)+p_angle(2)+p_angle(3))
        EllipsoidAngleX = 0;
        EllipsoidAngleY = EllipsoidAngleYmin + (EllipsoidAngleYmax-EllipsoidAngleYmin).*rand(1,1);
        EllipsoidAngleZ = 0;
    else
        EllipsoidAngleX = EllipsoidAngleXmin + (EllipsoidAngleXmax-EllipsoidAngleXmin).*rand(1,1);
        EllipsoidAngleY = EllipsoidAngleYmin + (EllipsoidAngleYmax-EllipsoidAngleYmin).*rand(1,1);
        EllipsoidAngleZ = 0;
    end
    EllipsoidAngleX = EllipsoidAngleX/180*pi;
    EllipsoidAngleY = EllipsoidAngleY/180*pi;
    EllipsoidAngleZ = EllipsoidAngleZ/180*pi;
    angles = [EllipsoidAngleX,EllipsoidAngleY,EllipsoidAngleZ];
       
    %椭球网格坐标生成
    [x,y,z] = ellipsoid(EllipsoidX,EllipsoidY,EllipsoidZ,EllipsoidAxisX,EllipsoidAxisY,EllipsoidAxisZ,Ellipsoid_precision);
    
    %旋转矩阵
    [RotationX,RotationY,RotationZ] = RotationMatrix(EllipsoidAngleX,EllipsoidAngleY,EllipsoidAngleZ);
                         
    %旋转操作
    for i = 1:size(x,1)
        for j = 1:size(x,2)
            xyz = RotationX*RotationY*RotationZ*[x(i,j) y(i,j) z(i,j) 1]';
            x(i,j) = xyz(1);
            y(i,j) = xyz(2);
            z(i,j) = xyz(3);
        end
    end
    [F,V] = surf2patch( x, y, z,'triangles');        %非常重要，输出stl的基础
    
    %扩大椭球的计算
    %计算随机扩大的范围
    EllipsoidDistanceX = EllipsoidDistanceXmin + (EllipsoidDistanceXmax-EllipsoidDistanceXmin)*rand(1,1);
    EllipsoidDistanceY = EllipsoidDistanceYmin + (EllipsoidDistanceYmax-EllipsoidDistanceYmin)*rand(1,1);
    EllipsoidDistanceZ = EllipsoidDistanceZmin + (EllipsoidDistanceZmax-EllipsoidDistanceZmin)*rand(1,1);
    
    %扩大之后的半轴长度
    EllipsoidAxisXExtend = EllipsoidAxisX + EllipsoidDistanceX/2;
    EllipsoidAxisYExtend = EllipsoidAxisY + EllipsoidDistanceY/2;
    EllipsoidAxisZExtend = EllipsoidAxisZ + EllipsoidDistanceZ/2;
    
    [x1,y1,z1] = ellipsoid(EllipsoidX,EllipsoidY,EllipsoidZ,EllipsoidAxisXExtend,EllipsoidAxisYExtend,EllipsoidAxisZExtend,Ellipsoid_precision);
    for i = 1:size(x1,1)
        for j = 1:size(x1,2)
            xyz = RotationX*RotationY*RotationZ*[x1(i,j) y1(i,j) z1(i,j) 1]';
            x1(i,j) = xyz(1);
            y1(i,j) = xyz(2);
            z1(i,j) = xyz(3);
        end
    end                                                                                         %计算扩大后椭球上的坐标

    %计算扩大后椭球的判别矩阵
    A = [ 1/EllipsoidAxisXExtend^2             0                          0               0
                     0               1/EllipsoidAxisYExtend^2             0               0
                     0                         0               1/EllipsoidAxisZExtend^2   0
                     0                         0                          0               -1 ];
                 
    T1 = [0 0 0 -EllipsoidX
          0 0 0 -EllipsoidY
          0 0 0 -EllipsoidZ
          0 0 0       1     ];
      
    [RotationX1,RotationY1,RotationZ1] = RotationMatrix(-EllipsoidAngleX,-EllipsoidAngleY,-EllipsoidAngleZ);
    
    B = T1*RotationZ1*RotationY1*RotationX1;
    
    P = B'*A*B;
        
    PXYZ = [x1(:),y1(:),z1(:),ones(numel(x1),1)]; 