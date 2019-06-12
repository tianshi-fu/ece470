function A = InvK3(T_1in0, M, S)
rng('shuffle')
[~,joints] = size(S);

switch joints
    case 1
    B1 = adjo(M^-1)*S(:,1);        
    theta = randn(joints,1);
    TC = M*expm(ssV(B1,theta(1)));
    ssVb = logm((TC^-1)*T_1in0);
    Vb = [ssVb(7); -ssVb(3); ssVb(2); ssVb(1:3,4)];
    n = norm(Vb);
    while n > 0.01
        n = norm(Vb);
        B1 = adjo(M^-1)*S(:,1);          
        J1 = B1;
        J = J1;
        mu = 0.1;
        thetadot = inv(J' * J + mu * eye(joints)) * J' * Vb;
        theta = theta + thetadot;
        TC = M*expm(ssV(B1,theta(1)));
        ssVb = logm((TC^-1)*T_1in0);
        Vb = [ssVb(7); -ssVb(3); ssVb(2); ssVb(1:3,4)];
    end
    case 2
        B1 = adjo(M^-1)*S(:,1);
        B2 = adjo(M^-1)*S(:,2);        
        theta = randn(joints,1);
        TC = M*expm(ssV(B1,theta(1)))*expm(ssV(B2,theta(2)));
        ssVb = logm((TC^-1)*T_1in0);
        Vb = [ssVb(7); -ssVb(3); ssVb(2); ssVb(1:3,4)];
        n = norm(Vb);
        while n > 0.01
            n = norm(Vb);
            B1 = adjo(M^-1)*S(:,1);
            B2 = adjo(M^-1)*S(:,2);           
            J2 = B2;
            J1 = adjo(expm(-ssV(B2,theta(2))))*B1;           
            J = [J1 J2];
            mu = 0.1;
            thetadot = inv(J' * J + mu * eye(joints)) * J' * Vb;
            theta = theta + thetadot;
            TC = M*expm(ssV(B1,theta(1)))*expm(ssV(B2,theta(2)));
            ssVb = logm((TC^-1)*T_1in0);
            Vb = [ssVb(7); -ssVb(3); ssVb(2); ssVb(1:3,4)];
        end
    case 3
        B1 = adjo(M^-1)*S(:,1);
        B2 = adjo(M^-1)*S(:,2);
        B3 = adjo(M^-1)*S(:,3);
        theta = randn(joints,1);
        TC = M*expm(ssV(B1,theta(1)))*expm(ssV(B2,theta(2)))*expm(ssV(B3,theta(3)));
        ssVb = logm((TC^-1)*T_1in0);
        Vb = [ssVb(7); -ssVb(3); ssVb(2); ssVb(1:3,4)];
        n = norm(Vb);
        while n > 0.01
            n = norm(Vb);
            B1 = adjo(M^-1)*S(:,1);
            B2 = adjo(M^-1)*S(:,2);
            B3 = adjo(M^-1)*S(:,3);
            J3 = B3;
            J2 = adjo(expm(-ssV(B3,theta(3))))*B2;
            J1 = adjo(expm(-ssV(B3,theta(3)))*expm(-ssV(B2,theta(2))))*B1;
            J = [J1 J2 J3];
            mu = 0.1;
            thetadot = inv(J' * J + mu * eye(joints)) * J' * Vb;
            theta = theta + thetadot;
            TC = M*expm(ssV(B1,theta(1)))*expm(ssV(B2,theta(2)))*expm(ssV(B3,theta(3)));
            ssVb = logm((TC^-1)*T_1in0);
            Vb = [ssVb(7); -ssVb(3); ssVb(2); ssVb(1:3,4)];
        end
    case 4
        B1 = adjo(M^-1)*S(:,1);
        B2 = adjo(M^-1)*S(:,2);
        B3 = adjo(M^-1)*S(:,3);
        B4 = adjo(M^-1)*S(:,4);
        theta = randn(joints,1);
        TC = M*expm(ssV(B1,theta(1)))*expm(ssV(B2,theta(2)))*expm(ssV(B3,theta(3)))*expm(ssV(B4,theta(4)));
        ssVb = logm((TC^-1)*T_1in0);
        Vb = [ssVb(7); -ssVb(3); ssVb(2); ssVb(1:3,4)];
        n = norm(Vb);
        while n > 0.01
            n = norm(Vb);
            B1 = adjo(M^-1)*S(:,1);
            B2 = adjo(M^-1)*S(:,2);
            B3 = adjo(M^-1)*S(:,3);
            B4 = adjo(M^-1)*S(:,4);
            J4 = B4;
            J3 = adjo(expm(-ssV(B4,theta(4))))*B3;
            J2 = adjo(expm(-ssV(B4,theta(4)))*expm(-ssV(B3,theta(3))))*B2;
            J1 = adjo(expm(-ssV(B4,theta(4)))*expm(-ssV(B3,theta(3)))*expm(-ssV(B2,theta(2))))*B1;
            J = [J1 J2 J3 J4];
            mu = 0.1;
            thetadot = inv(J' * J + mu * eye(joints)) * J' * Vb;
            theta = theta + thetadot;
            TC = M*expm(ssV(B1,theta(1)))*expm(ssV(B2,theta(2)))*expm(ssV(B3,theta(3)))*expm(ssV(B4,theta(4)));
            ssVb = logm((TC^-1)*T_1in0);
            Vb = [ssVb(7); -ssVb(3); ssVb(2); ssVb(1:3,4)];
        end
    case 5
        B1 = adjo(M^-1)*S(:,1);
        B2 = adjo(M^-1)*S(:,2);
        B3 = adjo(M^-1)*S(:,3);
        B4 = adjo(M^-1)*S(:,4);
        B5 = adjo(M^-1)*S(:,5);
        theta = randn(joints,1);
        TC = M*expm(ssV(B1,theta(1)))*expm(ssV(B2,theta(2)))*expm(ssV(B3,theta(3)))*expm(ssV(B4,theta(4)))*expm(ssV(B5,theta(5)));
        ssVb = logm((TC^-1)*T_1in0);
        Vb = [ssVb(7); -ssVb(3); ssVb(2); ssVb(1:3,4)];
        n = norm(Vb);
        while n > 0.01
            n = norm(Vb);
            B1 = adjo(M^-1)*S(:,1);
            B2 = adjo(M^-1)*S(:,2);
            B3 = adjo(M^-1)*S(:,3);
            B4 = adjo(M^-1)*S(:,4);
            B5 = adjo(M^-1)*S(:,5);
            J5 = B5;
            J4 = adjo(expm(-ssV(B5,theta(5))))*B4;
            J3 = adjo(expm(-ssV(B5,theta(5)))*expm(-ssV(B4,theta(4))))*B3;
            J2 = adjo(expm(-ssV(B5,theta(5)))*expm(-ssV(B4,theta(4)))*expm(-ssV(B3,theta(3))))*B2;
            J1 = adjo(expm(-ssV(B5,theta(5)))*expm(-ssV(B4,theta(4)))*expm(-ssV(B3,theta(3)))*expm(-ssV(B2,theta(2))))*B1;
            J = [J1 J2 J3 J4 J5];
            mu = 0.1;
            thetadot = inv(J' * J + mu * eye(joints)) * J' * Vb;
            theta = theta + thetadot;
            TC = M*expm(ssV(B1,theta(1)))*expm(ssV(B2,theta(2)))*expm(ssV(B3,theta(3)))*expm(ssV(B4,theta(4)))*expm(ssV(B5,theta(5)));
            ssVb = logm((TC^-1)*T_1in0);
            Vb = [ssVb(7); -ssVb(3); ssVb(2); ssVb(1:3,4)];
        end
    case 6
        B1 = adjo(M^-1)*S(:,1);
        B2 = adjo(M^-1)*S(:,2);
        B3 = adjo(M^-1)*S(:,3);
        B4 = adjo(M^-1)*S(:,4);
        B5 = adjo(M^-1)*S(:,5);
        B6 = adjo(M^-1)*S(:,6);
        theta = randn(joints,1);
        TC = M*expm(ssV(B1,theta(1)))*expm(ssV(B2,theta(2)))*expm(ssV(B3,theta(3)))*expm(ssV(B4,theta(4)))*expm(ssV(B5,theta(5)))*expm(ssV(B6,theta(6)));
        ssVb = logm((TC^-1)*T_1in0);
        Vb = [ssVb(7); -ssVb(3); ssVb(2); ssVb(1:3,4)];
        n = norm(Vb);
        while n > 0.01
            n = norm(Vb);
            B1 = adjo(M^-1)*S(:,1);
            B2 = adjo(M^-1)*S(:,2);
            B3 = adjo(M^-1)*S(:,3);
            B4 = adjo(M^-1)*S(:,4);
            B5 = adjo(M^-1)*S(:,5);
            B6 = adjo(M^-1)*S(:,6);
            J6 = B6;
            J5 = adjo(expm(-ssV(B6,theta(6))))*B5;
            J4 = adjo(expm(-ssV(B6,theta(6)))*expm(-ssV(B5,theta(5))))*B4;
            J3 = adjo(expm(-ssV(B6,theta(6)))*expm(-ssV(B5,theta(5)))*expm(-ssV(B4,theta(4))))*B3;
            J2 = adjo(expm(-ssV(B6,theta(6)))*expm(-ssV(B5,theta(5)))*expm(-ssV(B4,theta(4)))*expm(-ssV(B3,theta(3))))*B2;
            J1 = adjo(expm(-ssV(B6,theta(6)))*expm(-ssV(B5,theta(5)))*expm(-ssV(B4,theta(4)))*expm(-ssV(B3,theta(3)))*expm(-ssV(B2,theta(2))))*B1;
            J = [J1 J2 J3 J4 J5 J6];              
            mu = 0.1;
            thetadot = inv(J' * J + mu * eye(joints)) * J' * Vb;
            theta = theta + thetadot;
            TC = M*expm(ssV(B1,theta(1)))*expm(ssV(B2,theta(2)))*expm(ssV(B3,theta(3)))*expm(ssV(B4,theta(4)))*expm(ssV(B5,theta(5)))*expm(ssV(B6,theta(6)));
            ssVb = logm((TC^-1)*T_1in0);
            Vb = [ssVb(7); -ssVb(3); ssVb(2); ssVb(1:3,4)];
        end
    case 7
        B1 = adjo(M^-1)*S(:,1);
        B2 = adjo(M^-1)*S(:,2);
        B3 = adjo(M^-1)*S(:,3);
        B4 = adjo(M^-1)*S(:,4);
        B5 = adjo(M^-1)*S(:,5);
        B6 = adjo(M^-1)*S(:,6);
        B7 = adjo(M^-1)*S(:,7);
        theta = randn(joints,1);
        TC = M*expm(ssV(B1,theta(1)))*expm(ssV(B2,theta(2)))*expm(ssV(B3,theta(3)))*expm(ssV(B4,theta(4)))*expm(ssV(B5,theta(5)))*expm(ssV(B6,theta(6)))*expm(ssV(B7,theta(7)));
        ssVb = logm((TC^-1)*T_1in0);
        Vb = [ssVb(7); -ssVb(3); ssVb(2); ssVb(1:3,4)];
        n = norm(Vb);
        while n > 0.01
            n = norm(Vb);
            B1 = adjo(M^-1)*S(:,1);
            B2 = adjo(M^-1)*S(:,2);
            B3 = adjo(M^-1)*S(:,3);
            B4 = adjo(M^-1)*S(:,4);
            B5 = adjo(M^-1)*S(:,5);
            B6 = adjo(M^-1)*S(:,6);
            B7 = adjo(M^-1)*S(:,7);
            J7 = B7;
            J6 = adjo(expm(-ssV(B7,theta(7))))*B6;
            J5 = adjo(expm(-ssV(B7,theta(7)))*expm(-ssV(B6,theta(6))))*B5;
            J4 = adjo(expm(-ssV(B7,theta(7)))*expm(-ssV(B6,theta(6)))*expm(-ssV(B5,theta(5))))*B4;
            J3 = adjo(expm(-ssV(B7,theta(7)))*expm(-ssV(B6,theta(6)))*expm(-ssV(B5,theta(5)))*expm(-ssV(B4,theta(4))))*B3;
            J2 = adjo(expm(-ssV(B7,theta(7)))*expm(-ssV(B6,theta(6)))*expm(-ssV(B5,theta(5)))*expm(-ssV(B4,theta(4)))*expm(-ssV(B3,theta(3))))*B2;
            J1 = adjo(expm(-ssV(B7,theta(7)))*expm(-ssV(B6,theta(6)))*expm(-ssV(B5,theta(5)))*expm(-ssV(B4,theta(4)))*expm(-ssV(B3,theta(3)))*expm(-ssV(B2,theta(2))))*B1;
            J = [J1 J2 J3 J4 J5 J6 J7];              
            mu = 0.1;
            thetadot = inv(J' * J + mu * eye(joints)) * J' * Vb;
            theta = theta + thetadot;
            TC = M*expm(ssV(B1,theta(1)))*expm(ssV(B2,theta(2)))*expm(ssV(B3,theta(3)))*expm(ssV(B4,theta(4)))*expm(ssV(B5,theta(5)))*expm(ssV(B6,theta(6)))*expm(ssV(B7,theta(7)));
            ssVb = logm((TC^-1)*T_1in0);
            Vb = [ssVb(7); -ssVb(3); ssVb(2); ssVb(1:3,4)];
        end
    case 8
        B1 = adjo(M^-1)*S(:,1);
        B2 = adjo(M^-1)*S(:,2);
        B3 = adjo(M^-1)*S(:,3);
        B4 = adjo(M^-1)*S(:,4);
        B5 = adjo(M^-1)*S(:,5);
        B6 = adjo(M^-1)*S(:,6);
        B7 = adjo(M^-1)*S(:,7);
        B8 = adjo(M^-1)*S(:,8);
        theta = randn(joints,1);
        TC = M*expm(ssV(B1,theta(1)))*expm(ssV(B2,theta(2)))*expm(ssV(B3,theta(3)))*expm(ssV(B4,theta(4)))*expm(ssV(B5,theta(5)))*expm(ssV(B6,theta(6)))*expm(ssV(B7,theta(7)))*expm(ssV(B8,theta(8)));
        ssVb = logm((TC^-1)*T_1in0);
        Vb = [ssVb(7); -ssVb(3); ssVb(2); ssVb(1:3,4)];
        n = norm(Vb);
        while n > 0.01
            n = norm(Vb);
            B1 = adjo(M^-1)*S(:,1);
            B2 = adjo(M^-1)*S(:,2);
            B3 = adjo(M^-1)*S(:,3);
            B4 = adjo(M^-1)*S(:,4);
            B5 = adjo(M^-1)*S(:,5);
            B6 = adjo(M^-1)*S(:,6);
            B7 = adjo(M^-1)*S(:,7);
            B8 = adjo(M^-1)*S(:,8);
            J8 = B8;
            J7 = adjo(expm(-ssV(B8,theta(8))))*B7;
            J6 = adjo(expm(-ssV(B8,theta(8)))*expm(-ssV(B7,theta(7))))*B6;
            J5 = adjo(expm(-ssV(B8,theta(8)))*expm(-ssV(B7,theta(7)))*expm(-ssV(B6,theta(6))))*B5;
            J4 = adjo(expm(-ssV(B8,theta(8)))*expm(-ssV(B7,theta(7)))*expm(-ssV(B6,theta(6)))*expm(-ssV(B5,theta(5))))*B4;
            J3 = adjo(expm(-ssV(B8,theta(8)))*expm(-ssV(B7,theta(7)))*expm(-ssV(B6,theta(6)))*expm(-ssV(B5,theta(5)))*expm(-ssV(B4,theta(4))))*B3;
            J2 = adjo(expm(-ssV(B8,theta(8)))*expm(-ssV(B7,theta(7)))*expm(-ssV(B6,theta(6)))*expm(-ssV(B5,theta(5)))*expm(-ssV(B4,theta(4)))*expm(-ssV(B3,theta(3))))*B2;
            J1 = adjo(expm(-ssV(B8,theta(8)))*expm(-ssV(B7,theta(7)))*expm(-ssV(B6,theta(6)))*expm(-ssV(B5,theta(5)))*expm(-ssV(B4,theta(4)))*expm(-ssV(B3,theta(3)))*expm(-ssV(B2,theta(2))))*B1;
            J = [J1 J2 J3 J4 J5 J6 J7 J8];              
            mu = 0.1;
            thetadot = inv(J' * J + mu * eye(joints)) * J' * Vb;
            theta = theta + thetadot;
            TC = M*expm(ssV(B1,theta(1)))*expm(ssV(B2,theta(2)))*expm(ssV(B3,theta(3)))*expm(ssV(B4,theta(4)))*expm(ssV(B5,theta(5)))*expm(ssV(B6,theta(6)))*expm(ssV(B7,theta(7)))*expm(ssV(B8,theta(8)));
            ssVb = logm((TC^-1)*T_1in0);
            Vb = [ssVb(7); -ssVb(3); ssVb(2); ssVb(1:3,4)];
        end
    case 9
        B1 = adjo(M^-1)*S(:,1);
        B2 = adjo(M^-1)*S(:,2);
        B3 = adjo(M^-1)*S(:,3);
        B4 = adjo(M^-1)*S(:,4);
        B5 = adjo(M^-1)*S(:,5);
        B6 = adjo(M^-1)*S(:,6);
        B7 = adjo(M^-1)*S(:,7);
        B8 = adjo(M^-1)*S(:,8);
        B9 = adjo(M^-1)*S(:,9);
        theta = randn(joints,1);
        TC = M*expm(ssV(B1,theta(1)))*expm(ssV(B2,theta(2)))*expm(ssV(B3,theta(3)))*expm(ssV(B4,theta(4)))*expm(ssV(B5,theta(5)))*expm(ssV(B6,theta(6)))*expm(ssV(B7,theta(7)))*expm(ssV(B8,theta(8)))*expm(ssV(B9,theta(9)));
        ssVb = logm((TC^-1)*T_1in0);
        Vb = [ssVb(7); -ssVb(3); ssVb(2); ssVb(1:3,4)];
        n = norm(Vb);
        while n > 0.01
            n = norm(Vb);
            B1 = adjo(M^-1)*S(:,1);
            B2 = adjo(M^-1)*S(:,2);
            B3 = adjo(M^-1)*S(:,3);
            B4 = adjo(M^-1)*S(:,4);
            B5 = adjo(M^-1)*S(:,5);
            B6 = adjo(M^-1)*S(:,6);
            B7 = adjo(M^-1)*S(:,7);
            B8 = adjo(M^-1)*S(:,8);
            B9 = adjo(M^-1)*S(:,9);
            J9 = B9;
            J8 = adjo(expm(-ssV(B9,theta(9))))*B8;
            J7 = adjo(expm(-ssV(B9,theta(9)))*expm(-ssV(B8,theta(8))))*B7;
            J6 = adjo(expm(-ssV(B9,theta(9)))*expm(-ssV(B8,theta(8)))*expm(-ssV(B7,theta(7))))*B6;
            J5 = adjo(expm(-ssV(B9,theta(9)))*expm(-ssV(B8,theta(8)))*expm(-ssV(B7,theta(7)))*expm(-ssV(B6,theta(6))))*B5;
            J4 = adjo(expm(-ssV(B9,theta(9)))*expm(-ssV(B8,theta(8)))*expm(-ssV(B7,theta(7)))*expm(-ssV(B6,theta(6)))*expm(-ssV(B5,theta(5))))*B4;
            J3 = adjo(expm(-ssV(B9,theta(9)))*expm(-ssV(B8,theta(8)))*expm(-ssV(B7,theta(7)))*expm(-ssV(B6,theta(6)))*expm(-ssV(B5,theta(5)))*expm(-ssV(B4,theta(4))))*B3;
            J2 = adjo(expm(-ssV(B9,theta(9)))*expm(-ssV(B8,theta(8)))*expm(-ssV(B7,theta(7)))*expm(-ssV(B6,theta(6)))*expm(-ssV(B5,theta(5)))*expm(-ssV(B4,theta(4)))*expm(-ssV(B3,theta(3))))*B2;
            J1 = adjo(expm(-ssV(B9,theta(9)))*expm(-ssV(B8,theta(8)))*expm(-ssV(B7,theta(7)))*expm(-ssV(B6,theta(6)))*expm(-ssV(B5,theta(5)))*expm(-ssV(B4,theta(4)))*expm(-ssV(B3,theta(3)))*expm(-ssV(B2,theta(2))))*B1;
            J = [J1 J2 J3 J4 J5 J6 J7 J8 J9];              
            mu = 0.1;
            thetadot = inv(J' * J + mu * eye(joints)) * J' * Vb;
            theta = theta + thetadot;
            TC = M*expm(ssV(B1,theta(1)))*expm(ssV(B2,theta(2)))*expm(ssV(B3,theta(3)))*expm(ssV(B4,theta(4)))*expm(ssV(B5,theta(5)))*expm(ssV(B6,theta(6)))*expm(ssV(B7,theta(7)))*expm(ssV(B8,theta(8)))*expm(ssV(B9,theta(9)));
            ssVb = logm((TC^-1)*T_1in0);
            Vb = [ssVb(7); -ssVb(3); ssVb(2); ssVb(1:3,4)];
        end
    case 10
        B1 = adjo(M^-1)*S(:,1);
        B2 = adjo(M^-1)*S(:,2);
        B3 = adjo(M^-1)*S(:,3);
        B4 = adjo(M^-1)*S(:,4);
        B5 = adjo(M^-1)*S(:,5);
        B6 = adjo(M^-1)*S(:,6);
        B7 = adjo(M^-1)*S(:,7);
        B8 = adjo(M^-1)*S(:,8);
        B9 = adjo(M^-1)*S(:,9);
        B10 = adjo(M^-1)*S(:,10);
        theta = randn(joints,1);
        TC = M*expm(ssV(B1,theta(1)))*expm(ssV(B2,theta(2)))*expm(ssV(B3,theta(3)))*expm(ssV(B4,theta(4)))*expm(ssV(B5,theta(5)))*expm(ssV(B6,theta(6)))*expm(ssV(B7,theta(7)))*expm(ssV(B8,theta(8)))*expm(ssV(B9,theta(9)))*expm(ssV(B10,theta(10)));
        ssVb = logm((TC^-1)*T_1in0);
        Vb = [ssVb(7); -ssVb(3); ssVb(2); ssVb(1:3,4)];
        n = norm(Vb);
        while n > 0.01
            n = norm(Vb);
            B1 = adjo(M^-1)*S(:,1);
            B2 = adjo(M^-1)*S(:,2);
            B3 = adjo(M^-1)*S(:,3);
            B4 = adjo(M^-1)*S(:,4);
            B5 = adjo(M^-1)*S(:,5);
            B6 = adjo(M^-1)*S(:,6);
            B7 = adjo(M^-1)*S(:,7);
            B8 = adjo(M^-1)*S(:,8);
            B9 = adjo(M^-1)*S(:,9);
            B10 = adjo(M^-1)*S(:,10);
            J10 = B10;
            J9 = adjo(expm(-ssV(B10,theta(10))))*B9;
            J8 = adjo(expm(-ssV(B10,theta(10)))*expm(-ssV(B9,theta(9))))*B8;
            J7 = adjo(expm(-ssV(B10,theta(10)))*expm(-ssV(B9,theta(9)))*expm(-ssV(B8,theta(8))))*B7;
            J6 = adjo(expm(-ssV(B10,theta(10)))*expm(-ssV(B9,theta(9)))*expm(-ssV(B8,theta(8)))*expm(-ssV(B7,theta(7))))*B6;
            J5 = adjo(expm(-ssV(B10,theta(10)))*expm(-ssV(B9,theta(9)))*expm(-ssV(B8,theta(8)))*expm(-ssV(B7,theta(7)))*expm(-ssV(B6,theta(6))))*B5;
            J4 = adjo(expm(-ssV(B10,theta(10)))*expm(-ssV(B9,theta(9)))*expm(-ssV(B8,theta(8)))*expm(-ssV(B7,theta(7)))*expm(-ssV(B6,theta(6)))*expm(-ssV(B5,theta(5))))*B4;
            J3 = adjo(expm(-ssV(B10,theta(10)))*expm(-ssV(B9,theta(9)))*expm(-ssV(B8,theta(8)))*expm(-ssV(B7,theta(7)))*expm(-ssV(B6,theta(6)))*expm(-ssV(B5,theta(5)))*expm(-ssV(B4,theta(4))))*B3;
            J2 = adjo(expm(-ssV(B10,theta(10)))*expm(-ssV(B9,theta(9)))*expm(-ssV(B8,theta(8)))*expm(-ssV(B7,theta(7)))*expm(-ssV(B6,theta(6)))*expm(-ssV(B5,theta(5)))*expm(-ssV(B4,theta(4)))*expm(-ssV(B3,theta(3))))*B2;
            J1 = adjo(expm(-ssV(B10,theta(10)))*expm(-ssV(B9,theta(9)))*expm(-ssV(B8,theta(8)))*expm(-ssV(B7,theta(7)))*expm(-ssV(B6,theta(6)))*expm(-ssV(B5,theta(5)))*expm(-ssV(B4,theta(4)))*expm(-ssV(B3,theta(3)))*expm(-ssV(B2,theta(2))))*B1;
            J = [J1 J2 J3 J4 J5 J6 J7 J8 J9 J10];              
            mu = 0.1;
            thetadot = inv(J' * J + mu * eye(joints)) * J' * Vb;
            theta = theta + thetadot;
            TC = M*expm(ssV(B1,theta(1)))*expm(ssV(B2,theta(2)))*expm(ssV(B3,theta(3)))*expm(ssV(B4,theta(4)))*expm(ssV(B5,theta(5)))*expm(ssV(B6,theta(6)))*expm(ssV(B7,theta(7)))*expm(ssV(B8,theta(8)))*expm(ssV(B9,theta(9)))*expm(ssV(B10,theta(10)));
            ssVb = logm((TC^-1)*T_1in0);
            Vb = [ssVb(7); -ssVb(3); ssVb(2); ssVb(1:3,4)];
        end
    case 11
        B1 = adjo(M^-1)*S(:,1);
        B2 = adjo(M^-1)*S(:,2);
        B3 = adjo(M^-1)*S(:,3);
        B4 = adjo(M^-1)*S(:,4);
        B5 = adjo(M^-1)*S(:,5);
        B6 = adjo(M^-1)*S(:,6);
        B7 = adjo(M^-1)*S(:,7);
        B8 = adjo(M^-1)*S(:,8);
        B9 = adjo(M^-1)*S(:,9);
        B10 = adjo(M^-1)*S(:,10);
        B11 = adjo(M^-1)*S(:,11);
        theta = randn(joints,1);
        TC = M*expm(ssV(B1,theta(1)))*expm(ssV(B2,theta(2)))*expm(ssV(B3,theta(3)))*expm(ssV(B4,theta(4)))*expm(ssV(B5,theta(5)))*expm(ssV(B6,theta(6)))*expm(ssV(B7,theta(7)))*expm(ssV(B8,theta(8)))*expm(ssV(B9,theta(9)))*expm(ssV(B10,theta(10)))*expm(ssV(B11,theta(11)));
        ssVb = logm((TC^-1)*T_1in0);
        Vb = [ssVb(7); -ssVb(3); ssVb(2); ssVb(1:3,4)];
        n = norm(Vb);
        while n > 0.01
            n = norm(Vb);
            B1 = adjo(M^-1)*S(:,1);
            B2 = adjo(M^-1)*S(:,2);
            B3 = adjo(M^-1)*S(:,3);
            B4 = adjo(M^-1)*S(:,4);
            B5 = adjo(M^-1)*S(:,5);
            B6 = adjo(M^-1)*S(:,6);
            B7 = adjo(M^-1)*S(:,7);
            B8 = adjo(M^-1)*S(:,8);
            B9 = adjo(M^-1)*S(:,9);
            B10 = adjo(M^-1)*S(:,10);
            B11 = adjo(M^-1)*S(:,11);
            J11 = B11;
            J10 = adjo(expm(-ssV(B11,theta(11))))*B10;
            J9 = adjo(expm(-ssV(B11,theta(11)))*expm(-ssV(B10,theta(10))))*B9;
            J8 = adjo(expm(-ssV(B11,theta(11)))*expm(-ssV(B10,theta(10)))*expm(-ssV(B9,theta(9))))*B8;
            J7 = adjo(expm(-ssV(B11,theta(11)))*expm(-ssV(B10,theta(10)))*expm(-ssV(B9,theta(9)))*expm(-ssV(B8,theta(8))))*B7;
            J6 = adjo(expm(-ssV(B11,theta(11)))*expm(-ssV(B10,theta(10)))*expm(-ssV(B9,theta(9)))*expm(-ssV(B8,theta(8)))*expm(-ssV(B7,theta(7))))*B6;
            J5 = adjo(expm(-ssV(B11,theta(11)))*expm(-ssV(B10,theta(10)))*expm(-ssV(B9,theta(9)))*expm(-ssV(B8,theta(8)))*expm(-ssV(B7,theta(7)))*expm(-ssV(B6,theta(6))))*B5;
            J4 = adjo(expm(-ssV(B11,theta(11)))*expm(-ssV(B10,theta(10)))*expm(-ssV(B9,theta(9)))*expm(-ssV(B8,theta(8)))*expm(-ssV(B7,theta(7)))*expm(-ssV(B6,theta(6)))*expm(-ssV(B5,theta(5))))*B4;
            J3 = adjo(expm(-ssV(B11,theta(11)))*expm(-ssV(B10,theta(10)))*expm(-ssV(B9,theta(9)))*expm(-ssV(B8,theta(8)))*expm(-ssV(B7,theta(7)))*expm(-ssV(B6,theta(6)))*expm(-ssV(B5,theta(5)))*expm(-ssV(B4,theta(4))))*B3;
            J2 = adjo(expm(-ssV(B11,theta(11)))*expm(-ssV(B10,theta(10)))*expm(-ssV(B9,theta(9)))*expm(-ssV(B8,theta(8)))*expm(-ssV(B7,theta(7)))*expm(-ssV(B6,theta(6)))*expm(-ssV(B5,theta(5)))*expm(-ssV(B4,theta(4)))*expm(-ssV(B3,theta(3))))*B2;
            J1 = adjo(expm(-ssV(B11,theta(11)))*expm(-ssV(B10,theta(10)))*expm(-ssV(B9,theta(9)))*expm(-ssV(B8,theta(8)))*expm(-ssV(B7,theta(7)))*expm(-ssV(B6,theta(6)))*expm(-ssV(B5,theta(5)))*expm(-ssV(B4,theta(4)))*expm(-ssV(B3,theta(3)))*expm(-ssV(B2,theta(2))))*B1;
            J = [J1 J2 J3 J4 J5 J6 J7 J8 J9 J10 J11];              
            mu = 0.1;
            thetadot = inv(J' * J + mu * eye(joints)) * J' * Vb;
            theta = theta + thetadot;
            TC = M*expm(ssV(B1,theta(1)))*expm(ssV(B2,theta(2)))*expm(ssV(B3,theta(3)))*expm(ssV(B4,theta(4)))*expm(ssV(B5,theta(5)))*expm(ssV(B6,theta(6)))*expm(ssV(B7,theta(7)))*expm(ssV(B8,theta(8)))*expm(ssV(B9,theta(9)))*expm(ssV(B10,theta(10)))*expm(ssV(B11,theta(11)));
            ssVb = logm((TC^-1)*T_1in0);
            Vb = [ssVb(7); -ssVb(3); ssVb(2); ssVb(1:3,4)];
        end
    case 12
        B1 = adjo(M^-1)*S(:,1);
        B2 = adjo(M^-1)*S(:,2);
        B3 = adjo(M^-1)*S(:,3);
        B4 = adjo(M^-1)*S(:,4);
        B5 = adjo(M^-1)*S(:,5);
        B6 = adjo(M^-1)*S(:,6);
        B7 = adjo(M^-1)*S(:,7);
        B8 = adjo(M^-1)*S(:,8);
        B9 = adjo(M^-1)*S(:,9);
        B10 = adjo(M^-1)*S(:,10);
        B11 = adjo(M^-1)*S(:,11);
        B12 = adjo(M^-1)*S(:,12);
        theta = randn(joints,1);
        TC = M*expm(ssV(B1,theta(1)))*expm(ssV(B2,theta(2)))*expm(ssV(B3,theta(3)))*expm(ssV(B4,theta(4)))*expm(ssV(B5,theta(5)))*expm(ssV(B6,theta(6)))*expm(ssV(B7,theta(7)))*expm(ssV(B8,theta(8)))*expm(ssV(B9,theta(9)))*expm(ssV(B10,theta(10)))*expm(ssV(B11,theta(11)))*expm(ssV(B12,theta(12)));
        ssVb = logm((TC^-1)*T_1in0);
        Vb = [ssVb(7); -ssVb(3); ssVb(2); ssVb(1:3,4)];
        n = norm(Vb);
        while n > 0.01
            n = norm(Vb);
            B1 = adjo(M^-1)*S(:,1);
            B2 = adjo(M^-1)*S(:,2);
            B3 = adjo(M^-1)*S(:,3);
            B4 = adjo(M^-1)*S(:,4);
            B5 = adjo(M^-1)*S(:,5);
            B6 = adjo(M^-1)*S(:,6);
            B7 = adjo(M^-1)*S(:,7);
            B8 = adjo(M^-1)*S(:,8);
            B9 = adjo(M^-1)*S(:,9);
            B10 = adjo(M^-1)*S(:,10);
            B11 = adjo(M^-1)*S(:,11);
            B12 = adjo(M^-1)*S(:,12);
            J12 = B12;
            J11 = adjo(expm(-ssV(B12,theta(12))))*B11;
            J10 = adjo(expm(-ssV(B12,theta(12)))*expm(-ssV(B11,theta(11))))*B10;
            J9 = adjo(expm(-ssV(B12,theta(12)))*expm(-ssV(B11,theta(11)))*expm(-ssV(B10,theta(10))))*B9;
            J8 = adjo(expm(-ssV(B12,theta(12)))*expm(-ssV(B11,theta(11)))*expm(-ssV(B10,theta(10)))*expm(-ssV(B9,theta(9))))*B8;
            J7 = adjo(expm(-ssV(B12,theta(12)))*expm(-ssV(B11,theta(11)))*expm(-ssV(B10,theta(10)))*expm(-ssV(B9,theta(9)))*expm(-ssV(B8,theta(8))))*B7;
            J6 = adjo(expm(-ssV(B12,theta(12)))*expm(-ssV(B11,theta(11)))*expm(-ssV(B10,theta(10)))*expm(-ssV(B9,theta(9)))*expm(-ssV(B8,theta(8)))*expm(-ssV(B7,theta(7))))*B6;
            J5 = adjo(expm(-ssV(B12,theta(12)))*expm(-ssV(B11,theta(11)))*expm(-ssV(B10,theta(10)))*expm(-ssV(B9,theta(9)))*expm(-ssV(B8,theta(8)))*expm(-ssV(B7,theta(7)))*expm(-ssV(B6,theta(6))))*B5;
            J4 = adjo(expm(-ssV(B12,theta(12)))*expm(-ssV(B11,theta(11)))*expm(-ssV(B10,theta(10)))*expm(-ssV(B9,theta(9)))*expm(-ssV(B8,theta(8)))*expm(-ssV(B7,theta(7)))*expm(-ssV(B6,theta(6)))*expm(-ssV(B5,theta(5))))*B4;
            J3 = adjo(expm(-ssV(B12,theta(12)))*expm(-ssV(B11,theta(11)))*expm(-ssV(B10,theta(10)))*expm(-ssV(B9,theta(9)))*expm(-ssV(B8,theta(8)))*expm(-ssV(B7,theta(7)))*expm(-ssV(B6,theta(6)))*expm(-ssV(B5,theta(5)))*expm(-ssV(B4,theta(4))))*B3;
            J2 = adjo(expm(-ssV(B12,theta(12)))*expm(-ssV(B11,theta(11)))*expm(-ssV(B10,theta(10)))*expm(-ssV(B9,theta(9)))*expm(-ssV(B8,theta(8)))*expm(-ssV(B7,theta(7)))*expm(-ssV(B6,theta(6)))*expm(-ssV(B5,theta(5)))*expm(-ssV(B4,theta(4)))*expm(-ssV(B3,theta(3))))*B2;
            J1 = adjo(expm(-ssV(B12,theta(12)))*expm(-ssV(B11,theta(11)))*expm(-ssV(B10,theta(10)))*expm(-ssV(B9,theta(9)))*expm(-ssV(B8,theta(8)))*expm(-ssV(B7,theta(7)))*expm(-ssV(B6,theta(6)))*expm(-ssV(B5,theta(5)))*expm(-ssV(B4,theta(4)))*expm(-ssV(B3,theta(3)))*expm(-ssV(B2,theta(2))))*B1;
            J = [J1 J2 J3 J4 J5 J6 J7 J8 J9 J10 J11 J12];              
            mu = 0.1;
            thetadot = inv(J' * J + mu * eye(joints)) * J' * Vb;
            theta = theta + thetadot;
            TC = M*expm(ssV(B1,theta(1)))*expm(ssV(B2,theta(2)))*expm(ssV(B3,theta(3)))*expm(ssV(B4,theta(4)))*expm(ssV(B5,theta(5)))*expm(ssV(B6,theta(6)))*expm(ssV(B7,theta(7)))*expm(ssV(B8,theta(8)))*expm(ssV(B9,theta(9)))*expm(ssV(B10,theta(10)))*expm(ssV(B11,theta(11)))*expm(ssV(B12,theta(12)));
            ssVb = logm((TC^-1)*T_1in0);
            Vb = [ssVb(7); -ssVb(3); ssVb(2); ssVb(1:3,4)];
        end
    otherwise
        disp('Not implemented!!')
end




A = theta;
end