% intersection rayon Ur avec interface suivante (parabole)
function lambda1 = intersection_rayon_interface(a,b,c,vect_entering_point,vect_ray)
    
%     if (abs(a) < 1e-10)			% cas particulier d'une parabole qui est une droite
%         gamma = b * vect_entering_point(1) + c - vect_entering_point(2);
%         beta = b * vect_ray(1) - vect_ray(2);
%         lambda1 = -gamma / beta; % lambda1 is length of ray, the zero of linear function of lambda1 is searched
%     
%     else					% vraie parabole
%         gamma = a * vect_entering_point(1)^2 + b * vect_entering_point(1) + c - vect_entering_point(2);
%         if (abs(vect_ray(1)) < 1e-10) 
%             lambda1 = gamma / vect_ray(2);
%         
%         else
%             alpha = a * vect_ray(1)^2;
%             beta = 2 * a * vect_entering_point(1) * vect_ray(1) + b * vect_ray(1) - vect_ray(2);
%             delta = beta^2 - 4 * alpha * gamma;
%             delta = sqrt(delta);
%             lambda1 = (-beta - delta) / (2 * alpha);
%             lambda2 = (-beta + delta) / (2 * alpha);
%             if lambda1 <= 0
%                 lambda1 = lambda2; % lambda1 (or lambda2) is length of ray, the zeros of quadratic function of lambda1 are searched
%             elseif (lambda2<lambda1)&&(lambda2>0)
%                 lambda1 = lambda2; % lambda1 (or lambda2) is length of ray, the zeros of quadratic function of lambda1 are searched
%             end
%         end
%     end
%     
%     if isreal(lambda1)==0
%         lambda1 = 0;
%     end
    



	alpha = a * vect_ray(1)^2;
	beta = 2 * a * vect_entering_point(1) * vect_ray(1) + b * vect_ray(1) - vect_ray(2);
	gamma = a * vect_entering_point(1)^2 + b * vect_entering_point(1) + c - vect_entering_point(2);
    if (abs(alpha) < 1e-10) % straight line
		lambda1 = - gamma / beta;			% lambda1 is the length of the ray from Pt to interface
    else							% real parabola
		delta = beta^2 - 4. * alpha * gamma;
		delta = sqrt(delta);
		lambda1 = 2. * gamma / (- beta - delta);
		lambda2 = 2. * gamma / (- beta + delta);
        if lambda1 <= 0
            lambda1 = lambda2; % lambda1 (or lambda2) is length of ray, the zeros of quadratic function of lambda1 are searched
        elseif (lambda2<lambda1)&&(lambda2>0)
            lambda1 = lambda2; % lambda1 (or lambda2) is length of ray, the zeros of quadratic function of lambda1 are searched
        end
        
    end    
    if isreal(lambda1)==0
        lambda1 = 0;
    end
        