{
  description = "VSCode with default extensions and settings.";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixos-unstable";

    my-templates = {
      url = "github:manuelbb-upb/nix-templates";
      inputs.nixpkgs.follows = "nixpkgs";
    };

    nix-vscode-extensions = {
      url = "github:nix-community/nix-vscode-extensions";
      inputs.nixpkgs.follows = "nixpkgs";
    };
  };

  outputs = {
    self,
    nixpkgs,
    nix-vscode-extensions,
    my-templates,
  }@inputs:
  let
    system = "x86_64-linux";
    pkgs = nixpkgs.legacyPackages.${system};
    otherPkgs = my-templates.nix-matlab-ld.packages.${system};
    inherit (nix-vscode-extensions.extensions.${system}) vscode-marketplace;
    inherit (my-templates.packages.${system}) vscodium-with-default-extensions vscodium-local;

    winPkgs = pkgs.pkgsCross.mingwW64.buildPackages;

    additional-extensions-nixpkgs = with pkgs.vscode-extensions; [
      # extensions from nixpkgs
      ms-vscode.cpptools
      julialang.language-julia
      fortran-lang.linter-gfortran
    ];

    additional-extensions-marketplace = with vscode-marketplace; [
      # extensions from marketplace crawl
    ];

    additional-extensions-manual = pkgs.vscode-utils.extensionsFromVscodeMarketplace [
      # {
      #   name = "language-matlab";
      #   publisher = "MathWorks";
      #     version = "1.2.7";
      #       hash = "sha256-q/pjVGZPwKflZ0DSb95BUwePd4t5rCsXK+gyVFkJk7g=";
      # }
    ];

    vscode-base = pkgs.vscodium;
    
    vscode-with-exts = vscodium-with-default-extensions.override (prev: {
      vscode = vscode-base;
      vscodeExtensions = (
        prev.vscodeExtensions ++
        additional-extensions-nixpkgs ++ 
        additional-extensions-marketplace ++ 
        additional-extensions-manual
      );
    });

    vscode = vscodium-local.override {
      vscodium = vscode-with-exts;
      settings = {};  # right-merged with default-settings
      gitignore = "";
    };
  in
  {
    devShells.${system}.default = pkgs.mkShell {
      packages = (with pkgs; [
        cmake
        winPkgs.gcc
        winPkgs.gfortran
        gfortran
        fortls
      ]) ++ ( with otherPkgs; [
        (julia.override {
          julia-version = "1.11.1";
        })
      ]) ++[
        vscode
      ];
    };
  };
}
