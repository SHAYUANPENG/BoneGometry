%==========================================================================================
% Written by SAM,SDU��contact me if you have any question with this program   QQ��651587940    
%==========================================================================================
filename = 'BoneGeometry.stl'; 
%===============================================================================
%------------------------------------ ����ȷ�� ---------------------------------
%===============================================================================
BoneInnerRadius = 20;  % ��ͷ�ڰ뾶����ȷ����ֵ
BoneOuterRadius = 40;  % ��ͷ��뾶����ȷ����ֵ
BoneHeight = 100;      % ��ͷ�߶ȣ���ȷ����ֵ
precision = 40;        % Բ��һ���ж��ٸ��㣬���precision�������ģ�;��ȣ���ȷ����ֵ
radius_precision = 11; % �����������������ȷ����ֵ
Ellipsoid_precision = 40; %���򾫶ȣ�Խ�󣬾���Խ�ã���ȷ����ֵ

%����ʮ�㣬��Ϊ���ֺ�ȣ�����ĺ��֮��Ϊ���������޸Ĳ��Ͳ���

Thickness1 = 1.2;
Thickness2 = 0.8;
ThicknessSum = Thickness1 + Thickness2;
Cengshu = 10;

%�������ڲ��ڵĸ��ʣ�����Ҫ���޸�
p = 0.7;
%�����Ƿ���ת�ĸ��ʣ�����Ҫ���޸�
p_angle1 = 0.4;    %��������ת�ĸ�����0.4
p_angle2 = 0.2;    %ֻ��theta1
p_angle3 = 0.2;    %ֻ��theta2
p_angle4 = 0.2;    %����theta1����theta2
p_angle = [p_angle1,p_angle2,p_angle3,p_angle4];

%���±��������������min�Ĵ�����Сֵ����max�Ĵ������ֵ
%��X�ı�ʾX����Y��ʾY����Z��ʾZ����
%����X��Y��Z����������С��󳤶ȣ���ȷ�����������Χ
EllipsoidAxisXmin =  0.2;
EllipsoidAxisXmax =  0.4;
EllipsoidAxisYmin =  0.2;
EllipsoidAxisYmax =  0.4;
EllipsoidAxisZmin =  1;
EllipsoidAxisZmax =  3;  

%������X��Y��Z������ת�н���С��󳤶ȣ���ȷ�����������Χ
EllipsoidAngleXmin =  0;
EllipsoidAngleXmax =  5;
EllipsoidAngleYmin =  0;
EllipsoidAngleYmax =  5;
EllipsoidAngleZmin =  0;
EllipsoidAngleZmax =  0;      

%��������֮�����С���룬��ȷ�����������Χ
EllipsoidDistanceXmin = 0.1;
EllipsoidDistanceXmax = 0.5;
EllipsoidDistanceYmin = 0.1;
EllipsoidDistanceYmax = 0.5;
EllipsoidDistanceZmin = 3;
EllipsoidDistanceZmax = 15;  

%�������ĵ���ֵ�������ȷ�����������Χ
%���������Χ��Ҫ��Ϊ�˷�ֹ���������̫��Բ���߽��λ��
radius_min = BoneInnerRadius + 3*EllipsoidAxisXmax; %�뾶����
radius_max = BoneOuterRadius - 3*EllipsoidAxisXmax;
Height_min = 0 + 3*EllipsoidAxisZmax;               %�߶ȷ���
Height_max = BoneHeight - 3*EllipsoidAxisZmax;   

%������С������������ȷ�����������Χ
numberOfEllipsoidsMin = 10 ;
numberOfEllipsoidsMax = 30 ;


%======================================================================================
%-----------------------  ����ȫ�ֱ���������Ҫ���޸�  --------------------------------
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
%-------------------------------- Բ��������������� --------------------------------
%=====================================================================================

[Xinner,Yinner,Zinner] = cylinder(BoneInnerRadius, precision);           %��Բ�������꣬��ʱ�߶�Ϊ1
Zinner = Zinner*BoneHeight;                                              %�ڲ�Բ����߶�
[Finner,Vinner] = surf2patch(Xinner,Yinner,Zinner,'triangles');          %��ת��Ϊ�����˽ṹ��F��V������û�취��װ

[Xouter,Youter,Zouter] = cylinder(BoneOuterRadius, precision);           %��Բ�������꣬��ʱ�߶�Ϊ1
Zouter = Zouter*BoneHeight;                                              %���Բ����߶�
[Fouter,Vouter] = surf2patch(Xouter,Youter,Zouter,'triangles');

radius = linspace(BoneInnerRadius,BoneOuterRadius,radius_precision);     %������������
theta = (pi/180)*linspace(0,360,precision+1); % ����
[R,T] = meshgrid(radius,theta); % ����
XtopAndBottom = R.*cos(T); 
YtopAndBottom = R.*sin(T);      %ת��Ϊ�ѿ�������
Ztop = BoneHeight;              %����Z����
Zbottom = 0;                    %����Z����
[Ftop,Vtop] = surf2patch(XtopAndBottom, YtopAndBottom, BoneHeight.*ones(size(XtopAndBottom,1), size(XtopAndBottom,2)),'triangles');
[Fbottom,Vbottom] = surf2patch( XtopAndBottom, YtopAndBottom, 0.*XtopAndBottom,'triangles');


%======================================================================================
%----------------------------- �ڲ���������������� ----------------------------------      
%======================================================================================

% ���������������������������������̾����þ��ȷֲ�
numberOfEllipsoids = numberOfEllipsoidsMin + round((numberOfEllipsoidsMax-numberOfEllipsoidsMin).*rand(1,1)); 
Ellipsoids.numberOfEllipsoids = numberOfEllipsoids;

%��������
iEllipsoid = 1 ;
while iEllipsoid <= numberOfEllipsoids
    
    %��������
    [P,PXYZ,centroid,axisLength,angles,F,V,label] = CreatEllipsoid;

    %�б������Ƿ��ǵ�һ�����򣬷����ж������Ƿ������������ཻ
    if iEllipsoid == 1
        flag = zeros(1,size(PXYZ,1));
        for i = 1:size(PXYZ,1)    
             flag(i) = PXYZ(i,1)^2 + PXYZ(i,2)^2;
        end
        flag1 = flag(flag>BoneOuterRadius^2); flag2 =  flag(flag<BoneInnerRadius^2); 
        if isempty(flag1) && isempty(flag2)
           Ellipsoids.F(iEllipsoid).Faces  = F;
           Ellipsoids.V(iEllipsoid).vertex  = V;
           Ellipsoids.P(iEllipsoid).Matrix = P;  %ÿ��������б����
           Ellipsoids.PXYZ(iEllipsoid).coordinate = PXYZ;  %���������ϵĵ�
           Ellipsoids.centroid(iEllipsoid).coordinate = centroid; %�������ĵ�
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
           Ellipsoids.P(iEllipsoid).Matrix = P;  %ÿ��������б���󣬸��������������б�
           Ellipsoids.PXYZ(iEllipsoid).coordinate = PXYZ;  %���������ϵĵ�
           Ellipsoids.centroid(iEllipsoid).coordinate = centroid; %�������ĵ�
           Ellipsoids.axisLength(iEllipsoid).coordinate = axisLength;
           Ellipsoids.angles(iEllipsoid).coordinate = angles;  
           Ellipsoids.label(iEllipsoid) = label;
           iEllipsoid = iEllipsoid + 1;          
        end
    end
end
% ͳ�����ڲ�䣬���ڲ��ڵ���������
Tongji = Ellipsoids.label;
cengjianshu = sum(Tongji);
cengneishu = numberOfEllipsoids - cengjianshu;
fprintf('���������%d\n����������%d\n��������%d\n',cengjianshu,cengneishu,numberOfEllipsoids)
%=================================================================================
% -------------------------------- stl��ʽ���  ---------------------------------
%=================================================================================
% ��װF��V
% ����װԲ�������
Fzong = [Fbottom;Fouter+length(Vbottom(:,1))];
Vzong = [Vbottom;Vouter];
Fzong1 = [Fzong;Ftop+length(Vzong(:,1))];
Vzong1 = [Vzong;Vtop];
Fzong2 = [Fzong1;Finner+length(Vzong1(:,1))];
Vzong2 = [Vzong1;Vinner];

%ȥ���ظ��ڵ�
[Vzong2,indm,indn] = unique(Vzong2,'rows');           
Fzong2 = indn(Fzong2);

%��װ����F��V
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

% ��װԲ�����������
Fzong3 = [Fzong2; Fellipsoids+length(Vzong2(:,1))];
Vzong3 = [Vzong2; Vellipsoids];

%���stl
stlWrite(filename , Fzong3 , Vzong3);










    
    
    
    



