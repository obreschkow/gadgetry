#' @title Package environment
#'
#' @description Environment used to store default parameters used in \code{\link{writeparam}}, as well as optional other global variables. The default parameters are stored in \code{.gadgetry.env$default.parameters}.
#'
#' @seealso \code{\link{writeparam}}
#'
#' @author Danail Obreschkow
#'
#' @export

.gadgetry.env = new.env()

.gadgetry.env$default.parameters = list(OutputDir = '',
                                        SnapshotFileBase = 'snapshot',
                                        SnapFormat = 3,
                                        ICFormat = 1,
                                        InitCondFile = '',
                                        NumFilesPerSnapshot = 1,
                                        MaxFilesWithConcurrentIO = 0,
                                        TimeLimitCPU = 86400,
                                        CpuTimeBetRestartFile = 7200,
                                        MaxMemSize = 2000,
                                        TimeBegin = 0,
                                        TimeMax = 1,
                                        ComovingIntegrationOn = 0,
                                        BoxSize = 1,
                                        Omega0 = 0.3,
                                        OmegaLambda = 0.7,
                                        OmegaBaryon = 0.04,
                                        HubbleParam = 0.7,
                                        Hubble = 100,
                                        UnitLength_in_cm = 3.085678e21,
                                        UnitMass_in_g = 1.989e43,
                                        UnitVelocity_in_cm_per_s = 1e5,
                                        GravityConstantInternal = 0,
                                        TypeOfOpeningCriterion = 0,
                                        ErrTolTheta = 0.7,
                                        ErrTolForceAcc = 1,
                                        ErrTolThetaMax = 0.005,
                                        MaxSizeTimestep = 0.01,
                                        MinSizeTimestep = 0,
                                        ErrTolIntAccuracy = 0.025,
                                        CourantFac = 0.15,
                                        ActivePartFracForNewDomainDecomp = 0.01,
                                        TopNodeFactor = 2.5,
                                        OutputListOn = 0,
                                        OutputListFilename = 'output_times.txt',
                                        TimeOfFirstSnapshot = 0,
                                        TimeBetSnapshot = 0.1,
                                        TimeBetStatistics = 0.1,
                                        DesNumNgb = 64,
                                        MaxNumNgbDeviation = 2,
                                        InitGasTemp = 1e4,
                                        ArtBulkViscConst = 1,
                                        SofteningClassOfPartType0 = 0,
                                        SofteningClassOfPartType1 = 0,
                                        SofteningClassOfPartType2 = 0,
                                        SofteningClassOfPartType3 = 0,
                                        SofteningClassOfPartType4 = 0,
                                        SofteningClassOfPartType5 = 0,
                                        SofteningComovingClass0 = 0.5,
                                        SofteningMaxPhysClass0 = 0.5)
