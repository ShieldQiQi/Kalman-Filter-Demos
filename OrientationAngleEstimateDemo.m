% this demo is a simulation for using kalman filter in sensor fusion domain
% sensor: 3-axis gyroscope, 3-axis accelerometer and 3-axis magnetic sensor
% date: 2021-03-19
% Author: XiaoFangQi
%%
% define the 3-axis gyroscope sensor data
gryo_X = 0.0; gryo_Y = 0.0; gryo_Z = 0.1; 
% define the 3-axis accelerometer sensor data
acc_X = 5.0; acc_Y = 0.0; acc_Z = -8.4285; G=-9.8;
% define the 3-axis magnetic sensor data
magn_X = 4.0; magn_Y = 0.0; magn_Z = 0.0;
% 以下步骤计算磁场平行与地面的分量，该分量与重力方向垂直，可作为解算横摆角的依据
accVector = [acc_X; acc_Y; acc_Z];
magnVector = [magn_X; magn_Y; magn_Z];
magnVector = magnVector - dot(accVector,magnVector)*accVector/norm(accVector)/norm(accVector);
magn_X = magnVector(1); magn_Y = magnVector(2); magn_Z = magnVector(3); B = norm(magnVector);
%%
% simulate the orientation angles with reference to the EAST-NORTH-SKY coordinate
% get the initial orientation from the accelerometer & magnetic sensor
R_init = [1 0 0; 0 1 0; 0 0 1];
R_init(3,1) = acc_X/G; R_init(3,2) = acc_Y/G; R_init(3,3) = acc_Z/G;
R_init(2,1) = magn_X/B; R_init(2,2) = magn_Y/B; R_init(2,3) = magn_Z/B;
R_init(1,1:3) = cross(R_init(2,1:3),R_init(3,1:3));
% print the initial transformation
figure(1);
axis([-2 2 -2 2 -2 2]);
daspect([1 1 1]);
plotFrame([2, 0, 0; 0, 2 ,0; 0 ,0 ,2],'X','Y','Z','k','k','k');
hold on;
plotFrame(R_init,'x','y','z','r','g','b');

% update the orientation angle by integral of angular velocity
R_new = R_init;
% 设置积分间隔为0.01s
% 定义旋转矩阵微分角速度矩阵
pause(3);
% delta_t过大将导致积分误差增大，这个不仅仅是由速度不确定带来的误差，还有一部分为本身取极限不够极限带来的误差
delte_t = 0.001; 
OmegaVector = [gryo_X; gryo_Y; gryo_Z];
for t = 0:delte_t:15.7
    OmegaVectorInWorld = R_new*OmegaVector;
    % fprintf("%f %f %f \n",OmegaVectorInWorld(1),OmegaVectorInWorld(2),OmegaVectorInWorld(3));
    OmegaMatrix = [0 -OmegaVectorInWorld(3) OmegaVectorInWorld(2); OmegaVectorInWorld(3) 0 -OmegaVectorInWorld(1); -OmegaVectorInWorld(2) OmegaVectorInWorld(1) 0];
    R_new = (eye(3)+OmegaMatrix*delte_t)*R_new;
    if mod(t,delte_t*500) == 0 
        hold on;
        plotFrame(R_new,'x','y','z','r','g','b');
        pause(0.005);
    end
end
%%
function flag=plotFrame(Frame,lableX,lableY,lableZ,colorX,colorY,colorZ)
    title('kalman filter in sensor fusion');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    hold on;
    arrow3([0 0 0],[Frame(1,1) Frame(2,1) Frame(3,1)],colorX);
    hold on;
    text(Frame(1,1)-0.15,Frame(2,1)-0.15,Frame(3,1)-0.15,lableX);
    hold on;
    arrow3([0 0 0],[Frame(1,2) Frame(2,2) Frame(3,2)],colorY);
    hold on;
    text(Frame(1,2)-0.15,Frame(2,2)-0.15,Frame(3,2)-0.15,lableY);
    hold on;
    arrow3([0 0 0],[Frame(1,3) Frame(2,3) Frame(3,3)],colorZ);
    hold on;
    text(Frame(1,3)-0.15,Frame(2,3)-0.15,Frame(3,3)-0.15,lableZ);
    flag = 1;
end
