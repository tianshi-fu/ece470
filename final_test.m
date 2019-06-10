clear
clc
vrep=remApi('remoteApi');  %conection to VREP
clientID=vrep.simxStart('127.0.0.1',19997,true,true,5000,5);   %test the conection
if (clientID>-1)
 disp('Connected'); %code
    [returncode,targethandle]=vrep.simxGetObjectHandle(clientID,'ReferenceFrame',vrep.simx_opmode_blocking); %handles of relevant parts
    [returncode,joint1handle]=vrep.simxGetObjectHandle(clientID,'UR3_joint1',vrep.simx_opmode_blocking);
    [returncode,joint2handle]=vrep.simxGetObjectHandle(clientID,'UR3_joint2',vrep.simx_opmode_blocking);
    [returncode,joint3handle]=vrep.simxGetObjectHandle(clientID,'UR3_joint3',vrep.simx_opmode_blocking);
    [returncode,joint4handle]=vrep.simxGetObjectHandle(clientID,'UR3_joint4',vrep.simx_opmode_blocking);
    [returncode,joint5handle]=vrep.simxGetObjectHandle(clientID,'UR3_joint5',vrep.simx_opmode_blocking);
    [returncode,joint6handle]=vrep.simxGetObjectHandle(clientID,'UR3_joint6',vrep.simx_opmode_blocking);
    [returncode,EEhandle]=vrep.simxGetObjectHandle(clientID,'UR3_link7_visible',vrep.simx_opmode_blocking);
    [returncode,tiphandle]=vrep.simxGetObjectHandle(clientID,'Dummy_tip',vrep.simx_opmode_blocking);
    
    theta_init = [deg2rad(-90) deg2rad(-30) deg2rad(-45) 0 deg2rad(90) 0];
    old_theta = theta_init;
    
    vrep.simxSetJointTargetPosition(clientID,joint1handle,theta_init(1),vrep.simx_opmode_blocking); %set joint variables
    pause(0.05)
    vrep.simxSetJointTargetPosition(clientID,joint2handle,theta_init(2),vrep.simx_opmode_blocking);
    pause(0.05)
    vrep.simxSetJointTargetPosition(clientID,joint3handle,theta_init(3),vrep.simx_opmode_blocking);
    pause(0.05)
    vrep.simxSetJointTargetPosition(clientID,joint4handle,theta_init(4),vrep.simx_opmode_blocking);
    pause(0.05)
    vrep.simxSetJointTargetPosition(clientID,joint5handle,theta_init(5),vrep.simx_opmode_blocking);
    pause(0.05)
    vrep.simxSetJointTargetPosition(clientID,joint6handle,theta_init(6),vrep.simx_opmode_blocking);
    
    load('pathpoints_imported.mat')
    
    [returnCode,current_p]=vrep.simxGetObjectPosition(clientID,tiphandle,-1,vrep.simx_opmode_blocking);
    [returnCode,current_eulerAngles]=vrep.simxGetObjectOrientation(clientID,tiphandle,-1,vrep.simx_opmode_blocking);
    current_R = eul2rotm(current_eulerAngles,'XYZ');
    current_T = [current_R current_p'; 0 0 0 1];
    
    R_path = [0 -1 0; 0 0 -1; 1 0 0];   
    p_path = Pathpoints(1, 3:5)';
    next_T = [R_path p_path; 0 0 0 1];
    
    theta_c = IK_final(next_T,current_T,old_theta);
     
    vrep.simxSetJointTargetPosition(clientID,joint1handle,theta_init(1),vrep.simx_opmode_blocking); %set joint variables
    pause(0.05)
    vrep.simxSetJointTargetPosition(clientID,joint2handle,theta_init(2),vrep.simx_opmode_blocking);
    pause(0.05)
    vrep.simxSetJointTargetPosition(clientID,joint3handle,theta_init(3),vrep.simx_opmode_blocking);
    pause(0.05)
    vrep.simxSetJointTargetPosition(clientID,joint4handle,theta_init(4),vrep.simx_opmode_blocking);
    pause(0.05)
    vrep.simxSetJointTargetPosition(clientID,joint5handle,theta_init(5),vrep.simx_opmode_blocking);
    pause(0.05)
    vrep.simxSetJointTargetPosition(clientID,joint6handle,theta_init(6),vrep.simx_opmode_blocking);
    
    
    
else
    
    vrep.simxFinish(-1);  % finish the conection
end
vrep.delete(); % call the destructor!, finish the program
disp('Program ended');