FUNCTION INVERT_XTANHX(Y) RESULT (X0)
  ! Solve x * tanh(x) = y by dichotomy

  REAL, INTENT(IN) :: Y
  REAL             :: X0
  REAL             :: X_LOW, X_MEAN, X_UP
  REAL             :: EPS, STEP

  X_LOW = 0.0

  ! Find upper bound for x
  X_UP = 0.0
  STEP = MAX(Y, SQRT(Y))
  DO WHILE (X_UP*TANH(X_UP) < Y)
    X_UP = X_UP + STEP
  END DO

  ! Dichotomy
  EPS = 5.E-6
  DO WHILE (X_UP - X_LOW > EPS*X_UP)
    X_MEAN = (X_LOW + X_UP)/2
    IF (X_MEAN*TANH(X_MEAN) < Y) THEN
      X_LOW = X_MEAN
    ELSE
      X_UP = X_MEAN
    END IF
  END DO

  X0 = X_MEAN

  RETURN
END FUNCTION

