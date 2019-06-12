for f = 1:7
    n(f,f) = norm(p(:,f)-p(:,f+1));
    for g = 1:6
    n(f,f+1) = norm(p(:,f)-p(:,f+2));
    for h = 1:5
    n(f,f+2) = norm(p(:,f)-p(:,f+3));
    for k = 1:4
    n(f,f+3) = norm(p(:,f)-p(:,f+4));
    for l = 1:3
    n(f,f+4) = norm(p(:,f)-p(:,f+5));
    for m = 1:2
    n(f,f+5) = norm(p(:,f)-p(:,f+6));
    for o = 1
    n(f,f+6) = norm(p(:,f)-p(:,f+7));
    end
    end
    end
    end
    end
    end    
end
