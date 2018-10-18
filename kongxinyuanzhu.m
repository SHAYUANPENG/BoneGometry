%==========================================================================================
% Written by SAM,SDU，contact me if you have any question with this program   QQ：651587940    
%==========================================================================================
filename = 'BoneGeometry.stl'; 
%===============================================================================
%------------------------------------ 参数确定 ---------------------------------
%===============================================================================
BoneInnerRadius = 20;  % 骨头内半径，请确定数值
BoneOuterRadius = 40;  % 骨头外半径，请确定数值
BoneHeight = 100;      % 骨头高度，请确定数值
precision = 40;        % 圆柱一周有多少个点，提高precision可以提高模型精度，请确定数值
radius_precision = 11; % 径向网格点数量，请确定数值
Ellipsoid_precision = 40; %椭球精度，越大，精度越好，请确定数值

%共有十层，分为两种厚度，两层的厚度之和为常数，请修改层厚和层数

Thickness1 = 1.2;
Thickness2 = 0.8;
ThicknessSum = Thickness1 + Thickness2;
Cengshu = 10;

%椭球落在层内的概率，如需要请修改
p = 0.7;
%椭球是否旋转的概率，如需要请修改
p_angle1 = 0.4;    %不发生旋转的概率是0.4
p_angle2 = 0.2;    %只有theta1
p_angle3 = 0.2;    %只有theta2
p_angle4 = 0.2;    %既有theta1又有theta2
p_angle = [p_angle1,p_angle2,p_angle3,p_angle4];

%以下变量命名中最后是min的代表最小值，是max的代表最大值
%带X的表示X方向，Y表示Y方向，Z表示Z方向
%椭球X、Y、Z三个半轴最小最大长度，请确定随机变量范围
EllipsoidAxisXmin =  0.2;
EllipsoidAxisXmax =  0.4;
EllipsoidAxisYmin =  0.2;
EllipsoidAxisYmax =  0.4;
EllipsoidAxisZmin =  1;
EllipsoidAxisZmax =  3;  

%椭球绕X、Y、Z三轴旋转夹角最小最大长度，请确定随机变量范围
EllipsoidAngleXmin =  0;
EllipsoidAngleXmax =  5;
EllipsoidAngleYmin =  0;
EllipsoidAngleYmax =  5;
EllipsoidAngleZmin =  0;
EllipsoidAngleZmax =  0;      

%两个椭球之间的最小距离，请确定随机变量范围
EllipsoidDistanceXmin = 0.1;
EllipsoidDistanceXmax = 0.5;
EllipsoidDistanceYmin = 0.1;
EllipsoidDistanceYmax = 0.5;
EllipsoidDistanceZmin = 3;
EllipsoidDistanceZmax = 15;  

%椭球中心点出现的区域，请确定随机变量范围
%设置这个范围主要是为了防止椭球出现在太靠圆柱边界的位置
radius_min = BoneInnerRadius + 3*EllipsoidAxisXmax; %半径方向
radius_max = BoneOuterRadius - 3*EllipsoidAxisXmax;
Height_min = 0 + 3*EllipsoidAxisZmax;               %高度方向
Height_max = BoneHeight - 3*EllipsoidAxisZmax;   

%最大和最小椭球数量，请确定随机变量范围
numberOfEllipsoidsMin = 10 ;
numberOfEllipsoidsMax = 30 ;


%======================================================================================
%-----------------------  定义全局变量，无需要勿修改  --------------------------------
%======================================================================================
global Ellipsoids

Ellipsoids.BoneInnerRadius = BoneInnerRadius;
Ellipsoids.BoneOuterRadius = BoneOuterRadius;

Ellipsoids.Ellipsoid_precision = Ellipsoid_precision ;

Ellipsoids.Thickness1 = Thickness1;
Ellipsoids.Thickness2 = Thickness2;
Ellipsoids.Cengshu = Cengshu;
Ellipsoids.p = p;
Ellipsoids.p_angle.matrix = p_angle;

Ellipsoids.EllipsoidAxisXmin = EllipsoidAxisXmin;
Ellipsoids.EllipsoidAxisXmax = EllipsoidAxisXmax;
Ellipsoids.EllipsoidAxisYmin = EllipsoidAxisYmin;
Ellipsoids.EllipsoidAxisYmax = EllipsoidAxisYmax;
Ellipsoids.EllipsoidAxisZmin = EllipsoidAxisZmin;       
Ellipsoids.EllipsoidAxisZmax = EllipsoidAxisZmax;

Ellipsoids.EllipsoidAngleXmin = EllipsoidAngleXmin;
Ellipsoids.EllipsoidAngleXmax = EllipsoidAngleXmax;
Ellipsoids.EllipsoidAngleYmin = EllipsoidAngleYmin;
Ellipsoids.EllipsoidAngleYmax = EllipsoidAngleYmax;
Ellipsoids.EllipsoidAngleZmin = EllipsoidAngleZmin;
Ellipsoids.EllipsoidAngleZmax = EllipsoidAngleZmax;

Ellipsoids.EllipsoidDistanceXmin = EllipsoidDistanceXmin;
Ellipsoids.EllipsoidDistanceXmax = EllipsoidDistanceXmax;
Ellipsoids.EllipsoidDistanceYmin = EllipsoidDistanceYmin;
Ellipsoids.EllipsoidDistanceYmax = EllipsoidDistanceYmax;
Ellipsoids.EllipsoidDistanceZmin = EllipsoidDistanceZmin;
Ellipsoids.EllipsoidDistanceZmax = EllipsoidDistanceZmax;       

Ellipsoids.radius_min = radius_min;
Ellipsoids.radius_max = radius_max;
Ellipsoids.Height_min = Height_min;
Ellipsoids.Height_max = Height_max;   
%=====================================================================================
%-------------------------------- 圆柱外表面坐标生成 --------------------------------
%=====================================================================================

[Xinner,Yinner,Zinner] = cylinder(BoneInnerRadius, precision);           %内圆柱面坐标，此时高度为1
Zinner = Zinner*BoneHeight;                                              %内侧圆柱面高度
[Finner,Vinner] = surf2patch(Xinner,Yinner,Zinner,'triangles');          %不转化为有拓扑结构的F，V，后期没办法组装

[Xouter,Youter,Zouter] = cylinder(BoneOuterRadius, precision);           %外圆柱面坐标，此时高度为1
Zouter = Zouter*BoneHeight;                                              %外侧圆柱面高度
[Fouter,Vouter] = surf2patch(Xouter,Youter,Zouter,'triangles');

radius = linspace(BoneInnerRadius,BoneOuterRadius,radius_precision);     %径向网格数量
theta = (pi/180)*linspace(0,360,precision+1); % 周向
[R,T] = meshgrid(radius,theta); % 网格
XtopAndBottom = R.*cos(T); 
YtopAndBottom = R.*sin(T);      %转化为笛卡尔坐标
Ztop = BoneHeight;              %顶面Z坐标
Zbottom = 0;                    %底面Z坐标
[Ftop,Vtop] = surf2patch(XtopAndBottom, YtopAndBottom, BoneHeight.*ones(size(XtopAndBottom,1), size(XtopAndBottom,2)),'triangles');
[Fbottom,Vbottom] = surf2patch( XtopAndBottom, YtopAndBottom, 0.*XtopAndBottom,'triangles');


%======================================================================================
%----------------------------- 内部随机椭球坐标生成 ----------------------------------      
%======================================================================================

% 随机产生椭球数量，本程序所有随机过程均采用均匀分布
numberOfEllipsoids = numberOfEllipsoidsMin + round((numberOfEllipsoidsMax-numberOfEllipsoidsMin).*rand(1,1)); 
Ellipsoids.numberOfEllipsoids = numberOfEllipsoids;

%生成椭球
iEllipsoid = 1 ;
while iEllipsoid <= numberOfEllipsoids
    
    %创造椭球
    [P,PXYZ,centroid,axisLength,angles,F,V,label] = CreatEllipsoid;

    %判别椭球是否是第一个椭球，否则判断椭球是否与其他椭球相交
    if iEllipsoid == 1
        flag = zeros(1,size(PXYZ,1));
        for i = 1:size(PXYZ,1)    
             flag(i) = PXYZ(i,1)^2 + PXYZ(i,2)^2;
        end
        flag1 = flag(flag>BoneOuterRadius^2); flag2 =  flag(flag<BoneInnerRadius^2); 
        if isempty(flag1) && isempty(flag2)
           Ellipsoids.F(iEllipsoid).Faces  = F;
           Ellipsoids.V(iEllipsoid).vertex  = V;
           Ellipsoids.P(iEllipsoid).Matrix = P;  %每个椭球的判别矩阵
           Ellipsoids.PXYZ(iEllipsoid).coordinate = PXYZ;  %扩大椭球上的点
           Ellipsoids.centroid(iEllipsoid).coordinate = centroid; %椭球中心点
           Ellipsoids.axisLength(iEllipsoid).coordinate = axisLength;
           Ellipsoids.angles(iEllipsoid).coordinate = angles;    
           Ellipsoids.label(iEllipsoid) = label;
           iEllipsoid = iEllipsoid + 1;   
        end
    else
        flag = zeros(iEllipsoid-1,size(PXYZ,1));
        for isuccess = 1:iEllipsoid-1
             P1 = Ellipsoids.P(isuccess).Matrix;
             for i = 1:size(PXYZ,1)
                 flag(isuccess,i) = PXYZ(i,:)*P1*PXYZ(i,:)';  
             end
        end
        flag = flag(flag<0);
        
        flag1 = zeros(1,size(PXYZ,1));
        for i = 1:size(PXYZ,1)    
             flag1(i) = PXYZ(i,1)^2 + PXYZ(i,2)^2;
        end
        flag2 = flag1(flag1>BoneOuterRadius^2); flag3 =  flag1(flag1<BoneInnerRadius^2);
        if isempty(flag) && isempty(flag2) && isempty(flag3)
           Ellipsoids.F(iEllipsoid).Faces  = F;
           Ellipsoids.V(iEllipsoid).vertex  = V;
           Ellipsoids.P(iEllipsoid).Matrix = P;  %每个椭球的判别矩阵，根据扩大椭球来判别
           Ellipsoids.PXYZ(iEllipsoid).coordinate = PXYZ;  %扩大椭球上的点
           Ellipsoids.centroid(iEllipsoid).coordinate = centroid; %椭球中心点
           Ellipsoids.axisLength(iEllipsoid).coordinate = axisLength;
           Ellipsoids.angles(iEllipsoid).coordinate = angles;  
           Ellipsoids.label(iEllipsoid) = label;
           iEllipsoid = iEllipsoid + 1;          
        end
    end
end
% 统计落在层间，落在层内的椭球数量
Tongji = Ellipsoids.label;
cengjianshu = sum(Tongji);
cengneishu = numberOfEllipsoids - cengjianshu;
fprintf('层间椭球数%d\n层内椭球数%d\n椭球总数%d\n',cengjianshu,cengneishu,numberOfEllipsoids)
%=================================================================================
% -------------------------------- stl格式输出  ---------------------------------
%=================================================================================
% 组装F，V
% 先组装圆柱外表面
Fzong = [Fbottom;Fouter+length(Vbottom(:,1))];
Vzong = [Vbottom;Vouter];
Fzong1 = [Fzong;Ftop+length(Vzong(:,1))];
Vzong1 = [Vzong;Vtop];
Fzong2 = [Fzong1;Finner+length(Vzong1(:,1))];
Vzong2 = [Vzong1;Vinner];

%去掉重复节点
[Vzong2,indm,indn] = unique(Vzong2,'rows');           
Fzong2 = indn(Fzong2);

%组装椭球F，V
[x,y,z] = ellipsoid(0,0,0,1,2,3,Ellipsoid_precision);
[F,V] = surf2patch( x, y, z,'triangles');
constant = length(V(:,1));
mF = length(F(:,1));
Fellipsoids = zeros(mF*numberOfEllipsoids,3);
Vellipsoids = zeros(constant*numberOfEllipsoids,3);
for  iEllipsoid = 1: numberOfEllipsoids
     iFellipsoids = Ellipsoids.F(iEllipsoid).Faces + (iEllipsoid-1)*constant;
     Fellipsoids((1+(iEllipsoid-1)*mF):iEllipsoid*mF,:) = iFellipsoids;
     iVellipsoids = Ellipsoids.V(iEllipsoid).vertex;
     Vellipsoids((1+(iEllipsoid-1)*constant):iEllipsoid*constant,:) = iVellipsoids;
end

% 组装圆柱面和椭球面
Fzong3 = [Fzong2; Fellipsoids+length(Vzong2(:,1))];
Vzong3 = [Vzong2; Vellipsoids];

%输出stl
stlWrite(filename , Fzong3 , Vzong3);










    
    
    
    



